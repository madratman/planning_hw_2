/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include "mex.h"
#include <queue>

/* Input Arguments */
#define MAP_IN      prhs[0]
#define ARMSTART_IN prhs[1]
#define ARMGOAL_IN     prhs[2]
#define PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define PLAN_OUT    plhs[0]
#define PLANLENGTH_OUT  plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

// typedef for planar arm config - numofDOFs joint angles
typedef double* Config;

typedef struct {
    int X1, Y1;
    int X2, Y2;
    int Increment;
    int UsingYIndex;
    int DeltaX, DeltaY;
    int DTerm;
    int IncrE, IncrNE;
    int XIndex, YIndex;
    int Flipped;
} bresenham_param_t;

void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
    //take the nearest cell
    *pX = (int)(x/(double)(cellsize));
    if( x < 0) *pX = 0;
    if( *pX >= x_size) *pX = x_size-1;

    *pY = (int)(y/(double)(cellsize));
    if( y < 0) *pY = 0;
    if( *pY >= y_size) *pY = y_size-1;
}

void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
    params->UsingYIndex = 0;

    if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
        (params->UsingYIndex)++;

    if (params->UsingYIndex)
    {
        params->Y1=p1x;
        params->X1=p1y;
        params->Y2=p2x;
        params->X2=p2y;
    }
    else
    {
        params->X1=p1x;
        params->Y1=p1y;
        params->X2=p2x;
        params->Y2=p2y;
    }

     if ((p2x - p1x) * (p2y - p1y) < 0)
    {
        params->Flipped = 1;
        params->Y1 = -params->Y1;
        params->Y2 = -params->Y2;
    }
    else
        params->Flipped = 0;

    if (params->X2 > params->X1)
        params->Increment = 1;
    else
        params->Increment = -1;

    params->DeltaX=params->X2-params->X1;
    params->DeltaY=params->Y2-params->Y1;

    params->IncrE=2*params->DeltaY*params->Increment;
    params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
    params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

    params->XIndex = params->X1;
    params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
    if (params->UsingYIndex)
    {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
    else
    {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
    if (params->XIndex == params->X2)
    {
        return 0;
    }
    params->XIndex += params->Increment;
    if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
        params->DTerm += params->IncrE;
    else
    {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
    }
    return 1;
}

int IsValidLineSegment(double x0, double y0, double x1, double y1, double* map, int x_size, int y_size)

{
    bresenham_param_t params;
    int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
        
    //make sure the line segment is inside the environment
    if(x0 < 0 || x0 >= x_size ||
        x1 < 0 || x1 >= x_size ||
        y0 < 0 || y0 >= y_size ||
        y1 < 0 || y1 >= y_size)
        return 0;

    ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
    ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

    //iterate through the points on the segment
    get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
    do {
        get_current_point(&params, &nX, &nY);
        if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
                        return 0;
    } while (get_next_point(&params));

    return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double* map, int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
        
    //iterate through all the links starting with the base
    x1 = ((double)x_size)/2.0;
    y1 = 0;
    for(i = 0; i < numofDOFs; i++)
    {
        //compute the corresponding line segment
        x0 = x1;
        y0 = y1;
        x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
        y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

        //check the validity of the corresponding line segment
        if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
            return 0;
    }    
}

class Node
{
public:
    Node(Config config)
    {
        config_ = config;
        parent_idx_ = -1;
        numofDOFs_ = 5; // todo
        goal_thresh_joint_space_ = 10; // todo
        cost_ = 0;
    }

    Node(Config config, int parent_idx)
    {
        numofDOFs_ = 5; // todo
        config_ = config;
        parent_idx_ = parent_idx;
        cost_ = 0;
    }

    // constructor used for RRT*
    Node(Config config, int parent_idx, double cost)
    {
        numofDOFs_ = 5; // todo
        config_ = config;
        parent_idx_ = parent_idx;
        cost_ = cost;
    }

    void set_parent_idx(int parent_idx)
    {
        parent_idx_ = parent_idx;
    }

    double get_dist_to_config(Config config_2)
    {
        double total_dist = 0.0;
        double curr_diff = 0.0;
        double min_joint_angle = 0.0;
        for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
        {
            curr_diff = config_[joint_idx] - config_2[joint_idx];
            min_joint_angle = std::min(curr_diff, 2*M_PI - curr_diff);
            total_dist += pow(std::fabs(min_joint_angle), 2);
        }
        return sqrt(total_dist);
    }

    bool is_in_goal_region_joint_space(const Config &goal_config)
    {
        if (get_dist_to_config(goal_config) < goal_thresh_joint_space_)
        {
            return true;
        }
        return false;
    }

    // getters:
    int get_parent_idx()
    {
        return parent_idx_;
    }

    Config get_config()
    {
        return config_;
    }

    ~Node(){};

    // void set_numofDOFs(int numofDOFs)
    // {
    //     numofDOFs_ = numofDOFs;
    // }
// private:
    Config config_;// joint angles
    int parent_idx_;
    // static int const numofDOFs_;// = 4; initialize outside 
    // static double const goal_thresh_;
    int numofDOFs_;// = 4; initialize outside 
    double goal_thresh_joint_space_;
    double cost_; // for RRT* rewiring
};

class Tree
{
public:
    Tree(const Config &start_config,
            const Config &goal_config,
            double goal_thresh_cartesian,
            double epsilon,
            int num_steps_interp,
            int numofDOFs,
            double sample_goal_bias,
            double* map,
            int x_size,
            int y_size,
            double gamma, // rrt* user defined constant
            double epsilon_rrt_star,
            int tree_id)
    {
        start_config_ = start_config;
        goal_config_ = goal_config;
        num_samples_ = 0;
        is_goal_reached_ = false;
        goal_thresh_cartesian_ = goal_thresh_cartesian;
        epsilon_ = epsilon;
        numofDOFs_ = numofDOFs;
        sample_goal_bias_ = sample_goal_bias;
        num_steps_interp_ = num_steps_interp;
        map_ = map;
        map_x_size_ = x_size;
        map_y_size_ = y_size;
        goal_cartesian_ = new double[2];
        start_cartesian_ = new double[2];
        forward_kinematics(goal_config_, goal_cartesian_);
        forward_kinematics(start_config_, start_cartesian_);
        std::cout << "start config " << start_config_[0] << " " << start_config_[1]<< " " << start_config_[2]<< " " << start_config_[3] << std::endl;
        std::cout << "goal_config_ " << goal_config_[0] << " " << goal_config_[1]<< " " << goal_config_[2]<< " " << goal_config_[3] << std::endl;
        std::cout << "start_cartesian_ " << start_cartesian_[0] << ", " << start_cartesian_[1] << std::endl;
        std::cout << "goal_cartesian_ " << goal_cartesian_[0] << ", " << goal_cartesian_[1] << std::endl;
        tree_id_ = tree_id;
        insert_node(start_config);
        // RRT*
        delta_ = get_volume_R5_hyperball(); //todo generalize to numofDOFs_. currently hardcoded for R4 due to gamma function
        gamma_ = gamma;
        epsilon_rrt_star_ = epsilon_rrt_star;
    }

    ~Tree(){}; // todo free mem

    int get_nearest_index(const Config &config)
    {
        int idx_min = -1;
        double min_dist = std::numeric_limits<double>::infinity();
        for (int node_idx=0; node_idx < nodes_.size(); ++node_idx)
        {
            // std::cout << "node_idx : " << node_idx <<", ";
            // std::cout << "min_dist : " << min_dist <<", ";
            // std::cout << "current_dist : " << nodes_[node_idx].get_dist_to_config(config) << std::endl;
            if (nodes_[node_idx].get_dist_to_config(config) < min_dist)
            {
                min_dist = nodes_[node_idx].get_dist_to_config(config);
                idx_min = node_idx;
                // std::cout << "min_dist " << min_dist << ", idx_min " << idx_min << std::endl;
            }
        }
        // std::cout << std::endl;
        return idx_min;
    }

    // https://www.wikiwand.com/en/Volume_of_an_n-ball
    double get_volume_R5_hyperball()
    {
        return pow(M_PI, 2.5)/3.32335;//http://www.wolframalpha.com/input/?i=gamma+function(3.5)
    }

    // slide 36 @ http://www.cs.cmu.edu/~maxim/classes/robotplanning_grad/lectures/RRT_16782_fall17.pdf 
    void update_rrt_star_radius()
    {
        // double first_term = pow( gamma_ / delta_ * log(nodes_.size()) / nodes_.size() , 1.0/(double)numofDOFs_);
        if (nodes_.size() > 1)
            radius_rrt_star_ = std::min(pow( gamma_ / delta_ * log(nodes_.size()) / nodes_.size() , 1.0/(double)numofDOFs_), 
                                epsilon_rrt_star_);
            // radius_rrt_star_ = 5;
        return;
    }

    // nearest neighbours for RRT* 
    // TODO to make radius a param a tree, or not?
    void update_nearest_nodes_and_dist(Config sample_config)
    {
        nearest_nodes_dist_.clear();
        nearest_nodes_indices_.clear();
        closest_node_idx_ = 0;
        double min_dist = std::numeric_limits<double>::infinity();
        double curr_dist = 0.0;

        for (int node_idx; node_idx < nodes_.size(); node_idx++)
        {
            curr_dist = nodes_[node_idx].get_dist_to_config(sample_config);
            if (curr_dist < min_dist)
            {
                closest_node_idx_ = node_idx;
                min_dist = curr_dist;
                // std::cout << "min_dist " << min_dist  << "\n";
            }
            if (curr_dist <= radius_rrt_star_)
            {
                nearest_nodes_indices_.push_back(node_idx);
                nearest_nodes_dist_.push_back(curr_dist);
            }
        }
        closest_node_dist_ = min_dist;
        return;
    }

    // checks if transition is valid. 
    // TODO update RRT and RRT connect to use this
    // NOTE it checks till epsilon_ * sample_config !!!
    bool is_transition_valid(Config curr_config, Config sample_config)
    {
        Config intermediate_config = new double[numofDOFs_];
        for (int idx_step; idx_step < num_steps_interp_; idx_step++)
        {
            for (int joint_idx=0; joint_idx < numofDOFs_; joint_idx++)
            {
                // std::cout << "intermediate_config[ " << joint_idx << "] "<<  intermediate_config[joint_idx] 
                //             << " sample_config " << sample_config[joint_idx] << std::endl;
                intermediate_config[joint_idx] = curr_config[joint_idx] +
                            ((double)idx_step*(double)epsilon_*(sample_config[joint_idx] - curr_config[joint_idx])/(double)num_steps_interp_);
            }

            if (!IsValidArmConfiguration(intermediate_config, numofDOFs_, map_, map_x_size_, map_y_size_))
            { return false; }
        }
        return true;
    }

    bool extend_tree_rrt_star(Config sample_config)
    {
        std::cout << std::endl; 
        if (nodes_.size()==0)
        {
            Node start_node(start_config_, -1);
            nodes_.push_back(start_node);
            std::cout << "inserted first node" << std::endl;
            return true;
        }

        // steer sample to nearest neighbour and ensure it's within epsilon_rrt_star_ distance. 
        // TODO note this is implemented in a different way than RRT. 
        // I am not sure now which one is correct now/


        int idx_nearest = get_nearest_index(sample_config);
        Config intermediate_config = new double[numofDOFs_];
        for (int idx_step; idx_step < num_steps_interp_; idx_step++)
        {
            for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
            {
                intermediate_config[joint_idx] = nodes_[idx_nearest].get_config()[joint_idx] +
                            ((double)idx_step*(double)epsilon_*(sample_config[joint_idx] - nodes_[idx_nearest].get_config()[joint_idx])/(double)num_steps_interp_); 
            }

            if (!IsValidArmConfiguration(intermediate_config, numofDOFs_, map_, map_x_size_, map_y_size_))
            { return false; }                       
        }

        sample_config = intermediate_config;
        // idx_nearest = get_nearest_index(sample_config);
        update_rrt_star_radius();
        update_nearest_nodes_and_dist(sample_config);
        // closest_node_idx_ = idx_nearest;
        // std::cout << "num_nodes " << nodes_.size() << std::endl; 
        // std::cout << "radius_rrt_star_ " << radius_rrt_star_ << std::endl; 
        // std::cout << "closest_node_idx_ " << closest_node_idx_ << std::endl; 
        // std::cout << "idx_nearest " << closest_node_idx_ << std::endl; 

        // if(closest_node_dist_ > epsilon_rrt_star_)
        // {
        //     for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
        //     {
        //         // std::cout << "sample_config[" << joint_idx << "] " << sample_config[joint_idx] << std::endl;
        //         sample_config[joint_idx] = nodes_[closest_node_idx_].get_config()[joint_idx] +
        //                     epsilon_rrt_star_*((sample_config[joint_idx] - nodes_[closest_node_idx_].get_config()[joint_idx])/closest_node_dist_); 
        //         // std::cout << "sample_config[" << joint_idx << "] " << sample_config[joint_idx] << std::endl;
        //     }
        //     closest_node_dist_ = epsilon_rrt_star_;
        // }
        // Line 9 http://www.cs.cmu.edu/~maxim/classes/robotplanning_grad/lectures/RRT_16782_fall17.pdf
        if(is_transition_valid(nodes_[closest_node_idx_].get_config(), sample_config))
        {
            // dump collision checking results into vector to save computation for the second loop
            is_valid_nearest_nodes_vec_.clear();
            int min_node_idx = closest_node_idx_;
            double min_cost = nodes_[closest_node_idx_].cost_ + closest_node_dist_;
            double curr_cost = 0;
            bool is_transition_valid_curr = false;

            for (int loop_idx=0; loop_idx<nearest_nodes_indices_.size(); loop_idx++)
            {
                is_transition_valid_curr = is_transition_valid(nodes_[nearest_nodes_indices_[loop_idx]].get_config(),
                                                                sample_config);
                is_valid_nearest_nodes_vec_.push_back(is_transition_valid_curr);
                if (is_transition_valid_curr)
                {
                    // TODO bad unintuitive code. indices of indices, and indices together
                    curr_cost = nodes_[nearest_nodes_indices_[loop_idx]].cost_ + 
                                    nearest_nodes_dist_[loop_idx]; 
                    if (curr_cost < min_cost)
                    {
                        min_node_idx = nearest_nodes_indices_[loop_idx];
                        min_cost = curr_cost;
                    }
                }
            }

            // insert the sample into the tree with the currect parent idx and cost
            Node node_to_insert(sample_config, min_node_idx, min_cost);
            nodes_.push_back(node_to_insert);

            if (is_in_goal_region_cartesian_space(node_to_insert))
            {
                is_goal_reached_ = true;
                construct_path();
            }

            for (int loop_idx=0; loop_idx<nearest_nodes_indices_.size(); loop_idx++)
            {
                if (loop_idx == min_node_idx)
                    continue;
                curr_cost = nodes_[nearest_nodes_indices_[loop_idx]].cost_ 
                                + nearest_nodes_dist_[loop_idx]; 
                if (is_valid_nearest_nodes_vec_[loop_idx] && curr_cost < nearest_nodes_dist_[loop_idx])
                {
                    nodes_[nearest_nodes_indices_[loop_idx]].cost_ = curr_cost;
                    nodes_[nearest_nodes_indices_[loop_idx]].parent_idx_ = loop_idx;
                }
            }
        }
    }

    // todo confirm don't need y_size as we already check validity before adding to tree
    void forward_kinematics(const Config &config, double* cartesian)
    {
        double x0,y0,x1,y1;
        int i;
        //iterate through all the links starting with the base
        x1 = ((double)map_x_size_)/2.0;
        y1 = 0;
        for(i = 0; i < numofDOFs_; i++)
        {
            //compute the corresponding line segment
            x0 = x1;
            y0 = y1;
            x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-config[i]);
            y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-config[i]);
        }

        cartesian[0] = x1;
        cartesian[1] = y1;
    }

    bool is_in_goal_region_cartesian_space(Node& node)
    {
        double* node_cartesian = new double[numofDOFs_];
        forward_kinematics(node.config_, node_cartesian);
        double dist = pow(node_cartesian[0]-goal_cartesian_[0], 2) +
                      pow(node_cartesian[1]-goal_cartesian_[1], 2);
        dist = sqrt(dist);

        if (dist < goal_thresh_cartesian_)
        {   
            return true;
        }
        return false;
    }

    bool are_close_enough_cartesian_space(Config &config_1, Config &config_2)
    {
        double* cartesian_1 = new double[numofDOFs_];
        double* cartesian_2 = new double[numofDOFs_];
        forward_kinematics(config_1, cartesian_1);
        forward_kinematics(config_2, cartesian_2);
        std::cout << "cartesian_1 : (" << cartesian_1[0] << ", " << cartesian_1[1] << std::endl;
        std::cout << "cartesian_2 : (" << cartesian_2[0] << ", " << cartesian_2[1] << std::endl;
        double dist = pow(cartesian_1[0]-cartesian_2[0], 2) +
                      pow(cartesian_1[1]-cartesian_2[1], 2);
        dist = sqrt(dist);

        std::cout << "cartesian dist "<< dist << std::endl;
        if (dist < goal_thresh_cartesian_)
        {   
            return true;
        }
        return false;
    }


    bool insert_node(Config sample_config)
    {
        if (nodes_.size()==0)
        {
            Node start_node(start_config_, -1);
            nodes_.push_back(start_node);
            std::cout << "inserted first node" << std::endl;
            return true;
        }

        int idx_nearest = get_nearest_index(sample_config);
        // std::cout << "got nearest index " << idx_nearest << std::endl; 

        // interpolate to ensure path is valid
        Config intermediate_config = new double[numofDOFs_];
        for (int idx_step; idx_step < num_steps_interp_; idx_step++)
        {
            // std::cout<< "idx_step " << idx_step << std::endl;
            // interpolate the intermediate config
            for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
            {
                // std::cout << " nearest " << nodes_[idx_nearest].get_config()[joint_idx];
                // std::cout << ", sample " << sample_config[joint_idx];

                intermediate_config[joint_idx] = nodes_[idx_nearest].get_config()[joint_idx] +
                            ((double)idx_step*(double)epsilon_*(sample_config[joint_idx] - nodes_[idx_nearest].get_config()[joint_idx])/(double)num_steps_interp_); 
                // std::cout << ", to be added " << intermediate_config[joint_idx] << std::endl;
            }

            // add it to tree if it is valid 
            if (!IsValidArmConfiguration(intermediate_config, numofDOFs_, map_, map_x_size_, map_y_size_))
            { return false; }                       
        }

        Node new_node(intermediate_config, idx_nearest); // parent_idx is idx_nearest
        nodes_.push_back(new_node);
        idx_nearest = nodes_.size()-1; // nearest now points to last inserted
        // std::cout << "added node " <<std::endl;
        if (is_in_goal_region_cartesian_space(new_node))
        {
            is_goal_reached_ = true;
            construct_path();
            return true; // goal reached itself! => return true 
        }
        // managed to insert all intermediate configs => return true
        return true;
    }

    bool extend_tree(Config sample_config)
    {
        if (nodes_.size()==0)
        {
            Node start_node(start_config_, -1);
            nodes_.push_back(start_node);
            std::cout << "inserted first node" << std::endl;
            return true;
        }

        int idx_nearest = get_nearest_index(sample_config);
        // std::cout << "got nearest index " << idx_nearest << std::endl; 

        // interpolate to ensure path is valid
        Config intermediate_config = new double[numofDOFs_];
        for (int idx_step; idx_step < num_steps_interp_; idx_step++)
        {
            // std::cout<< "idx_step " << idx_step << std::endl;
            // interpolate the intermediate config
            for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
            {
                // std::cout << " nearest " << nodes_[idx_nearest].get_config()[joint_idx];
                // std::cout << ", sample " << sample_config[joint_idx];
                intermediate_config[joint_idx] = nodes_[idx_nearest].get_config()[joint_idx] +
                            ((double)idx_step*(double)epsilon_*(sample_config[joint_idx] - nodes_[idx_nearest].get_config()[joint_idx])/(double)num_steps_interp_); 
                // std::cout << ", to be added " << intermediate_config[joint_idx] << std::endl;
            }

            // add it to tree if it is valid 
            if (IsValidArmConfiguration(intermediate_config, numofDOFs_, map_, map_x_size_, map_y_size_))
            {
                // intermediate_config is now epsilon_ towards the sample_config
                Node new_node(intermediate_config, idx_nearest); // parent_idx is idx_nearest
                nodes_.push_back(new_node);
                idx_nearest = nodes_.size()-1; // nearest now points to last inserted
                // std::cout << "added node " <<std::endl;
                if (are_close_enough_cartesian_space(intermediate_config, sample_config))
                {            
                    construct_path();
                    return true;  
                }
            }
            else
            {
                // std::cout << "invalid config " << std::endl;
                // some intermediate config is in collision => return false
                return false;
            }                       
        }


        // managed to insert all intermediate configs => return true
        return true;
    }

    void construct_path()
    {
        path_.clear();
        int idx_waypt = nodes_.size()-1;
        while (idx_waypt!=-1)
        {
            path_.push_back(nodes_[idx_waypt].get_config());
            idx_waypt = nodes_[idx_waypt].get_parent_idx();
        }
        std::reverse(path_.begin(), path_.end());
    }

    Config get_random_sample()
    {
        // std::cout << "num_samples_ : " << num_samples_ << std::endl; 
        // std::cout << "num_nodes : " << nodes_.size() << std::endl; 
        num_samples_++;

        if ((double) rand()/(double) RAND_MAX < sample_goal_bias_)
        {
            // std::cout << "sampling goal config" << std::endl;
            return goal_config_;
        }
        else
        {
            Config sample_config = new double[numofDOFs_];

            for(int i = 0; i < numofDOFs_; ++i) 
            {
                sample_config[i] = (double)rand()*(2*M_PI)/(double) RAND_MAX;
                // std::cout << "sample_config[i]" << sample_config[i] << std::endl;
            }
            return sample_config;
        }
    } 

    std::vector<Config> get_path()
    {
        return path_;
    }

    Config start_config_;
    Config goal_config_;
    double* goal_cartesian_;
    double* start_cartesian_;
    std::vector<Node> nodes_;
    std::vector<Config> path_;
    int num_samples_;
    bool is_goal_reached_;
    int numofDOFs_;
    double* map_;
    int map_x_size_;
    int map_y_size_;
    double sample_goal_bias_;
    double epsilon_;
    double goal_thresh_cartesian_; // todo
    int num_steps_interp_;
    int tree_id_;
    // rrt* stuff:
    double delta_;
    double gamma_;
    double epsilon_rrt_star_;
    double radius_rrt_star_;
    int closest_node_idx_;
    std::vector<int> nearest_nodes_indices_;
    std::vector<double> nearest_nodes_dist_;
    std::vector<int> is_valid_nearest_nodes_vec_;
    double closest_node_dist_;

};

static void planner_rrt(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)


{
    *plan = NULL;
    *planlength = 0;
    // double epsilon = 0.01; // looks sexier
    double epsilon = 0.1; // runs faster
    int num_max_iters = (int)1e6;
    double sample_goal_bias = 0.3;
    double goal_thresh_cartesian = 5;
    int num_steps_interp = 20 ;
    int tree_id = 0 ;

    Tree tree(armstart_anglesV_rad, 
            armgoal_anglesV_rad, 
            goal_thresh_cartesian,
            epsilon,
            num_steps_interp,
            numofDOFs, 
            sample_goal_bias,
            map,
            x_size,
            y_size,
            0, // rrt* user defined constant
            0,
            tree_id);

    bool got_path=false;

    for (int iter=0; iter < num_max_iters; iter++)
    {
        if (tree.is_goal_reached_)
        {
            std::cout << " FOUND SOLUTION " << std::endl;
            got_path = true;
            break;
        }

        Config sample_config = tree.get_random_sample();
        if(!IsValidArmConfiguration(sample_config, numofDOFs, map, x_size, y_size))
            continue;
        tree.insert_node(sample_config);
    }

    if (got_path)
    {
        std::vector<Config> path_vector = tree.get_path();
        *planlength = (int)path_vector.size();
        *plan = (double**) malloc(path_vector.size()*sizeof(double*));

        for (int i = 0; i < path_vector.size(); i++)
        {
            (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
            for(int j = 0; j < numofDOFs; j++)
            {
                (*plan)[i][j] = path_vector[i][j];
            }
        }            
    }

    return;
}

static void planner_rrt_star(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    *plan = NULL;
    *planlength = 0;
    // double epsilon = 0.01; // looks sexier
    double epsilon = 0.4; // runs faster
    double sample_goal_bias = 0.5;
    double goal_thresh_cartesian = 5; // this is epsilon
    int num_steps_interp = 10;
    int tree_id = 0;
    //  user defined constant for RRT* radius
    // slide 36 @ http://www.cs.cmu.edu/~maxim/classes/robotplanning_grad/lectures/RRT_16782_fall17.pdf
    double gamma = 100;
    double epsilon_rrt_star = M_PI/4; 
    double radius; 
    int num_max_iters = 1e6;

    Tree tree(armstart_anglesV_rad, 
            armgoal_anglesV_rad, 
            goal_thresh_cartesian,
            epsilon,
            num_steps_interp,
            numofDOFs, 
            sample_goal_bias,
            map,
            x_size,
            y_size,
            0, // rrt* user defined constant
            0,
            tree_id);

    bool got_path=false;

    for (int iter=0; iter < num_max_iters; iter++)
    {
        // std::cout << "iter " << iter << std::endl;
        if (tree.is_goal_reached_)
        {
            std::cout << " FOUND SOLUTION " << std::endl;
            got_path = true;
            break;
        }

        Config sample_config = tree.get_random_sample();
        if(!IsValidArmConfiguration(sample_config, numofDOFs, map, x_size, y_size))
            continue;
        tree.extend_tree_rrt_star(sample_config);
    }

    if (got_path)
    {
        std::vector<Config> path_vector = tree.get_path();
        *planlength = (int)path_vector.size();
        *plan = (double**) malloc(path_vector.size()*sizeof(double*));

        for (int i = 0; i < path_vector.size(); i++)
        {
            (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
            for(int j = 0; j < numofDOFs; j++)
            {
                (*plan)[i][j] = path_vector[i][j];
            }
        }            
    }

    return;
}

// void swap_trees(Tree* current_tree_ptr, bool &is_start_tree, Tree &start_tree, Tree &goal_tree)
// {
//     if (is_start_tree)
//     {
//         current_tree_ptr = &goal_tree;
//         is_start_tree = false;
//     }
//     else
//     {
//         current_tree_ptr = &start_tree;
//         is_start_tree = true;
//     }
//     return;
// }

static void planner_rrt_connect(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    *plan = NULL;
    *planlength = 0;
    double epsilon = 0.1;
    double sample_goal_bias = 0.3;
    double goal_thresh_cartesian = 5;
    int num_steps_interp = 20;

    Tree start_tree(armstart_anglesV_rad, 
            armgoal_anglesV_rad, 
            goal_thresh_cartesian,
            epsilon,
            num_steps_interp,
            numofDOFs, 
            sample_goal_bias,
            map,
            x_size,
            y_size,
            0,
            0,
            0);


   Tree goal_tree(armgoal_anglesV_rad,
            armstart_anglesV_rad, 
            goal_thresh_cartesian,
            epsilon,
            num_steps_interp,
            numofDOFs, 
            sample_goal_bias,
            map,
            x_size,
            y_size,
            0,
            0,
            100);

    Tree* current_tree_ptr;
    current_tree_ptr = &start_tree;
    // current_tree_ptr = &goal_tree;
    bool is_start_tree = true;
    Config current_config; // todo unintuitive code. Config is already a pointer 
    Node* nearest_neighbour;
    bool sample_was_valid = false; 
    bool path_found = false; 
    bool trees_connect = false; 
    int num_max_iters = 1e6;

    for (int iter=0; iter < num_max_iters; iter++)
    {
        if (path_found)
        {
            break;
        }
        if (is_start_tree)
            std::cout << "Iter " << iter << ", Start Tree" << std::endl;
        else
            std::cout << "Iter " << iter << ", Goal Tree" << std::endl;
        std::cout << "start_tree size " << start_tree.nodes_.size() 
                  << ", goal_tree size " << goal_tree.nodes_.size() << std::endl; 

        // Get random sample
        Config sample_config = current_tree_ptr->get_random_sample();
        if(!IsValidArmConfiguration(sample_config, numofDOFs, map, x_size, y_size))
            continue;
        // try to insert in current tree
        sample_was_valid = current_tree_ptr->insert_node(sample_config);
        // std::cout << "sample_was_valid " << sample_was_valid << std::endl;

        if (current_tree_ptr->is_goal_reached_)
        {
            std::cout << "current_tree_ptr reached goal itself !!!" << std::endl;
            path_found = true;
        }
        // now we swap the trees. Note in the original paper, swap is written after the next
        // if clause, however we're trying to now "tree->insert()"" sample_config into the 
        // other tree for the "Connect" part of the algorithm
        // swap_trees(current_tree_ptr, is_start_tree, start_tree, goal_tree);
        if (is_start_tree)
        {
            current_tree_ptr = &goal_tree;
            is_start_tree = false;
        }
        else
        {
            current_tree_ptr = &start_tree;
            is_start_tree = true;
        }

        // if (sample_was_valid) is equivalent to Line 4, Figure 5 : if not(Extend(Ta, q_rand) = Trapped)
        // of the RRT Connect paper http://citeseerx.ist.psu.edu/viewdocdownload?doi=10.1.1.2.2971&rep=rep1&type=pdf
        if (sample_was_valid)
        {
            trees_connect = current_tree_ptr->extend_tree(sample_config);
            if ((start_tree.nodes_.size()>0) && (goal_tree.nodes_.size()>0))
            {
                std::cout << "checking if connect " << std::endl;
                if (trees_connect)
                {
                    std::cout << "trees connect !!!" << std::endl;
                    path_found = true;
                    break;
                }            
            }
        }
    }

    if (path_found)
    {
        start_tree.construct_path();
        goal_tree.construct_path();
        std::vector<Config> path_start_tree = start_tree.get_path();
        std::vector<Config> path_goal_tree = goal_tree.get_path(); 
        std::reverse(path_goal_tree.begin(), path_goal_tree.end());

        std::cout << "path_start_tree.size() " << path_start_tree.size() << std::endl;
        std::cout << "path_goal_tree.size() " << path_goal_tree.size() << std::endl;
        *planlength = (int)path_start_tree.size() + (int)path_goal_tree.size();
        *plan = (double**) malloc(*planlength*sizeof(double*));
        std::cout << "*planlength " << *planlength << std::endl;

        // fill in the start tree's path in plan
        if (path_start_tree.size() > 0)
        {
            for (int i = 0; i < path_start_tree.size(); i++)
            {
                (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
                for(int j = 0; j < numofDOFs; j++)
                {
                    (*plan)[i][j] = path_start_tree[i][j];
                }
            }            
        } 
        if (path_goal_tree.size() > 0)
        {
            for (int i = 0; i < path_goal_tree.size(); i++)
            {
                (*plan)[i+path_start_tree.size()] = (double*) malloc(numofDOFs*sizeof(double)); 
                for(int j = 0; j < numofDOFs; j++)
                {
                    (*plan)[i+path_start_tree.size()][j] = path_goal_tree[i][j];
                }
            }
        }
    }
    return;
}

class PRM_Node
{
public:
    PRM_Node(){};
    PRM_Node(Config config)
    {
        config_ = config;
        is_conn_to_start_ = false;
        is_conn_to_goal_ = false;
        idx_ = -1;
        numofDOFs_ = 5;
    }

    PRM_Node(Config config, 
            bool is_conn_to_start, 
            bool is_conn_to_goal)
    {
        config_ = config;
        is_conn_to_start_ = is_conn_to_start;
        is_conn_to_goal_ = is_conn_to_goal;
        idx_ = -1;
        numofDOFs_ = 5;
    }

    PRM_Node(Config config, 
            bool is_conn_to_start, 
            bool is_conn_to_goal,
            int idx)
    {
        config_ = config;
        is_conn_to_start_ = is_conn_to_start;
        is_conn_to_goal_ = is_conn_to_goal;
        idx_ = idx;
        numofDOFs_ = 5;
    }

    double get_dist_to_config(Config config_2)
    {
        double total_dist = 0.0;
        double curr_diff = 0.0;
        double min_joint_angle = 0.0;
        for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
        {
            curr_diff = config_[joint_idx] - config_2[joint_idx];
            min_joint_angle = std::min(curr_diff, 2*M_PI - curr_diff);
            total_dist += pow(std::fabs(min_joint_angle), 2);
        }
        return sqrt(total_dist);
    }

    ~PRM_Node(){};
    
    Config config_;
    // std::vector<PRM_Node> neighbour_vec_;
    std::vector<int> nearest_nodes_indices_;
    std::vector<double> nearest_nodes_dist_;
    bool is_conn_to_start_;
    bool is_conn_to_goal_;
    int idx_;
    int numofDOFs_;
    // int parent_idx_;
    PRM_Node* parent_ptr_;
};

class PRMGraph
{
public:
    PRMGraph(){};
    PRMGraph(Config start_config,
              Config goal_config,
              int numofDOFs,
              double epsilon,
              double epsilon_rrt_star,
              int num_steps_interp,
              double sample_goal_bias,
              int num_max_iters,
              double gamma,
              double* map,
              int x_size,
              int y_size)
    {
        PRM_Node start_node_(start_config, true, false, 1);
        PRM_Node goal_node_(goal_config, false, true, -1);
        nodes_.push_back(start_node_);
        nodes_.push_back(goal_node_);
        delta_ = get_volume_R5_hyperball();
        gamma_ = gamma;
        epsilon_ = epsilon;
        map_ = map;
        map_x_size_ = x_size;
        map_y_size_ = y_size;
        numofDOFs_ = numofDOFs;
        num_steps_interp_ = num_steps_interp;
        epsilon_rrt_star_ = epsilon_rrt_star;
        num_samples_ = 0;
        sample_goal_bias_ = sample_goal_bias;
        start_config_ = start_config;
        goal_config_ = goal_config;
        num_max_iters_= num_max_iters;
    }

    ~PRMGraph(){};

    double get_volume_R5_hyperball()
    {
        return pow(M_PI, 2.5)/3.32335;//http://www.wolframalpha.com/input/?i=gamma+function(3.5)
    }

    void update_rrt_star_radius()
    {
        if (nodes_.size() > 1)
            radius_rrt_star_ = std::min(pow( gamma_ / delta_ * log(nodes_.size()) / nodes_.size() , 1.0/(double)numofDOFs_), 
                                epsilon_rrt_star_);
        return;
    }

    Config get_random_sample()
    {
        if ((double) rand()/(double) RAND_MAX < sample_goal_bias_)
        {
            // std::cout << "sampling goal config" << std::endl;
            return goal_config_;
        }
        else
        {
            Config sample_config = new double[numofDOFs_];

            for(int i = 0; i < numofDOFs_; ++i) 
            {
                sample_config[i] = (double)rand()*(2*M_PI)/(double) RAND_MAX;
            }
            return sample_config;
        }
    }

    void update_nearest_nodes_and_dist(PRM_Node* sample_node_ptr)
    {
        // sample_node_ptr->nearest_nodes_dist.clear();
        // sample_node_ptr->nearest_nodes_indices.clear();
        int closest_node_idx = 0;
        double min_dist = std::numeric_limits<double>::infinity();
        double curr_dist = 0.0;

        for (int node_idx; node_idx < nodes_.size(); node_idx++)
        {
            curr_dist = nodes_[node_idx].get_dist_to_config(sample_node_ptr->config_);
            if (curr_dist < min_dist)
            {
                closest_node_idx = node_idx;
                min_dist = curr_dist;
            }
            if (curr_dist <= radius_rrt_star_)
            {
                sample_node_ptr->nearest_nodes_indices_.push_back(node_idx);
                sample_node_ptr->nearest_nodes_dist_.push_back(curr_dist);
            }
        }
        // closest_node_dist = min_dist;
        return;
    }

    bool is_transition_valid(Config curr_config, Config sample_config)
    {
        Config intermediate_config = new double[numofDOFs_];
        for (int idx_step; idx_step < num_steps_interp_; idx_step++)
        {
            for (int joint_idx=0; joint_idx < numofDOFs_; joint_idx++)
            {
                intermediate_config[joint_idx] = curr_config[joint_idx] +
                            ((double)idx_step*(double)epsilon_*(sample_config[joint_idx] - curr_config[joint_idx])/(double)num_steps_interp_);
            }

            if (!IsValidArmConfiguration(intermediate_config, numofDOFs_, map_, map_x_size_, map_y_size_))
            { return false; }
        }
        return true;
    }

    bool propagate_start_goal_connectivity(PRM_Node* curr_node_ptr) 
    {
        bool is_start_and_goal_conn = (curr_node_ptr->is_conn_to_start_
                                    && curr_node_ptr->is_conn_to_goal_);
        PRM_Node* curr_neighbour_ptr;
        bool is_curr_neighbour_conn_to_start;
        bool is_curr_neighbour_conn_to_goal;
        for (int loop_idx = 0; loop_idx < curr_node_ptr->nearest_nodes_indices_.size(); loop_idx++) 
        {
            double curr_neighbour_idx = curr_node_ptr->nearest_nodes_indices_[loop_idx];
            *curr_neighbour_ptr = nodes_[curr_neighbour_idx];
            bool curr_neighbour_was_conn_to_start = curr_neighbour_ptr->is_conn_to_start_;
            bool curr_neighbour_was_conn_to_goal = curr_neighbour_ptr->is_conn_to_goal_;
            curr_neighbour_ptr->is_conn_to_start_ = (curr_neighbour_ptr->is_conn_to_start_ 
                                                        || curr_node_ptr->is_conn_to_start_);

            curr_neighbour_ptr->is_conn_to_goal_ = (curr_neighbour_ptr->is_conn_to_goal_ 
                                                        || curr_node_ptr->is_conn_to_goal_);
            
            if (curr_neighbour_was_conn_to_start != curr_node_ptr->is_conn_to_start_ 
                    || curr_neighbour_was_conn_to_goal != curr_node_ptr->is_conn_to_goal_) 
            {
                if( propagate_start_goal_connectivity(curr_neighbour_ptr) ) 
                {
                    is_start_and_goal_conn = true;
                }
            }
        }
        return is_start_and_goal_conn;
    }

    void build_graph()
    {
        for (int iter=0; iter < num_max_iters_; iter++)
        {
            Config sample_config = get_random_sample();
            if(!IsValidArmConfiguration(sample_config, numofDOFs_, map_, map_x_size_, map_y_size_))
                continue;
            num_samples_++;
            update_rrt_star_radius();
            PRM_Node* curr_node_ptr = new PRM_Node(sample_config, false, false, -1);
            update_nearest_nodes_and_dist(curr_node_ptr);
            int curr_node_idx = nodes_.size()+1;
            nodes_.push_back(*curr_node_ptr);

            PRM_Node* curr_neighbour;
            for (int loop_idx=0; loop_idx < curr_node_ptr->nearest_nodes_indices_.size(); loop_idx++)
            {
                double curr_neighbour_idx = curr_node_ptr->nearest_nodes_indices_[loop_idx];
                if ( !is_transition_valid(nodes_[curr_neighbour_idx].config_, sample_config) );
                    { continue; }
                // add as neighbour of the curr neighbour!
                nodes_[curr_neighbour_idx].nearest_nodes_indices_.push_back(curr_node_idx);
                curr_node_ptr->nearest_nodes_indices_.push_back(curr_neighbour_idx);
            }

            if (propagate_start_goal_connectivity(curr_node_ptr))
                { break; }
        }
        return;
    }

    std::vector<PRM_Node> nodes_;
    PRM_Node start_node_;
    PRM_Node goal_node_;
    Config start_config_;
    Config goal_config_;
    int numofDOFs_;
    int num_steps_interp_;
    double epsilon_;
    double radius_rrt_star_;
    double epsilon_rrt_star_;
    double delta_;
    double gamma_;
    int map_x_size_;
    int map_y_size_;
    double* map_;
    int num_samples_;
    double sample_goal_bias_;
    int num_max_iters_;
};

static void planner_prm(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    *plan = NULL;
    *planlength = 0;
    // double epsilon = 0.01; // looks sexier
    double epsilon = 0.4; // runs faster
    double sample_goal_bias = 0.5;
    double goal_thresh_cartesian = 5; // this is epsilon
    int num_steps_interp = 10;
    int tree_id = 0;
    //  user defined constant for RRT* radius
    // slide 36 @ http://www.cs.cmu.edu/~maxim/classes/robotplanning_grad/lectures/RRT_16782_fall17.pdf
    double gamma = 100;
    double epsilon_rrt_star = M_PI/4; 
    double radius; 
    int num_max_iters = 100;

    bool got_path=false;
    PRMGraph prm_graph_obj;
    prm_graph_obj.PRMGraph(Config armstart_anglesV_rad,
                        Config armgoal_anglesV_rad,
                        int numofDOFs,
                        double epsilon,
                        double epsilon_rrt_star,
                        int num_steps_interp,
                        double sample_goal_bias,
                        int num_max_iters,
                        double gamma,
                        double* map,
                        int x_size,
                        int y_size);

    std::cout << prm_graph_obj.gamma_ << std::endl;
    // prm_graph_obj.build_graph();

    // std::queue<PRM_Node*> prm_queue;
    // prm_queue.push(prm_graph.start_node_);
            
    // PRM_Node* curr_node_ptr;
    // while(prm_queue.size() != 0) 
    // {
    //     curr_node_ptr = prm_queue.front();
    //     prm_queue.pop();
    //     if (*curr_node_ptr == prm_graph.goal_node_) 
    //     {
    //         std::cout << "SOLUTION FOUND!!!" << std::endl;
    //         break;
    //     }
    //     PRM_Node* curr_neighbour_ptr;
    //     for(int loop_idx = 0; loop_idx < curr_node_ptr->nearest_nodes_indices_.size(); loop_idx++)
    //     {
    //         double curr_neighbour_idx = curr_node_ptr->nearest_nodes_indices_[loop_idx];
    //         *curr_neighbour_ptr = nodes_[curr_neighbour_idx];
    //         if (curr_neighbour_ptr->idx == -1) 
    //         {
    //             curr_neighbour_ptr->parent_ptr_ = curr_node_ptr;
    //             curr_neighbour_ptr->idx = curr_node_ptr->idx + 1;
    //             prm_queue.push(curr_neighbour_ptr);
    //         }
    //     }
    // }

    // *plan = (double**) malloc(curr_node_ptr->idx * sizeof(double*));
    // *planlength = curr_node_ptr->idx;

    // for (int i = *planlength - 1; i >= 0; i--)
    // {
    //     (*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
    //     for(int j = 0; j < numofDOFs; j++)
    //     {
    //         (*plan)[i][j] = curr_node_ptr->config_[j];
    //     }
    //     curr_node_ptr = curr_node_ptr->parent_ptr_;
    // }
    

}
//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
            int nrhs, const mxArray*prhs[])
         
{ 
        
    /* Check for proper number of arguments */    
    if (nrhs != 4)
    { 
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                            "Four input arguments required."); 
    } 
    else if (nlhs != 2)
    {
        mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                            "One output argument required."); 
    } 
            
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1)
    {
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                            "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN)))
    {
                        mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                            "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);

    //get the planner id
    int planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if(planner_id < 0 || planner_id > 3)
    {
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidplanner_id",
                            "planner id should be between 0 and 3 inclusive");         
    }
    
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    
    //you can may be call the corresponding planner function here
    if (planner_id==RRT)
    {
        planner_rrt(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }
    if (planner_id==RRTCONNECT)
    {
        planner_rrt_connect(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }
    if (planner_id==RRTSTAR)
    {
        planner_rrt_star(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }
    
    //dummy planner which only computes interpolated path
    // planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    
    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
            plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*) mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    return;
}

