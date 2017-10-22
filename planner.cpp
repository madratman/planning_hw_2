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
    Node(const Config config)
    {
        config_ = config;
        parent_idx_ = -1;
        numofDOFs_ = 4; // todo
        goal_thresh_ = 10; // todo
    }

    Node(const Config config, int parent_idx)
    {
        config_ = config;
        parent_idx_ = parent_idx;
    }

    void set_parent_idx(int parent_idx)
    {
        parent_idx_ = parent_idx;
    }

    double get_dist_to_config(Config config_2)
    {
        double total_dist = 0.0;
        for (int joint_idx=0; joint_idx<numofDOFs_; joint_idx++)
        {
            total_dist = total_dist + 
                        std::abs(std::min(config_[joint_idx] - config_2[joint_idx],
                                2*M_PI-(config_[joint_idx] - config_2[joint_idx])));
        }
        return total_dist;
    }

    bool is_in_goal_region(const Config &goal_config)
    {
        if (get_dist_to_config(goal_config) < goal_thresh_)
        {
            return true;
        }
        else
        {
            return false;
        }
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

    void set_numofDOFs(int numofDOFs)
    {
        numofDOFs_ = numofDOFs;
    }
// private:
    Config config_;// joint angles
    int parent_idx_;
    // static int const numofDOFs_;// = 4; initialize outside 
    // static double const goal_thresh_;
    int numofDOFs_;// = 4; initialize outside 
    double goal_thresh_;
};


class Tree
{
public:
    Tree(Config start_config, 
            Config goal_config,
            int num_max_iters,
            double epsilon,
            int numofDOFs,
            double sample_goal_bias)
    {
        start_config_ = start_config;
        goal_config_ = goal_config;
        num_samples_ = 0;
        is_goal_reached_ = false;
        num_max_iters_ = num_max_iters;
        epsilon_ = epsilon;
        numofDOFs_ = numofDOFs;
        insert_node(start_config);
        sample_goal_bias_ = sample_goal_bias;
    }

    ~Tree(){}; // todo free mem

    int get_nearest_index(const Config &config)
    {
        int idx_min = -1;
        double min_dist = std::numeric_limits<double>::infinity();
        for (int node_idx=0; node_idx < nodes_.size(); node_idx++)
        {
            // std::cout << "node_idx : " << node_idx <<", ";
            // std::cout << "min_dist : " << min_dist <<", ";
            // std::cout << "current_dist : " << nodes_[node_idx].get_dist_to_config(config) << std::endl;
            if (nodes_[node_idx].get_dist_to_config(config) < min_dist)
            {
                // if (min_dist <= 0.01)
                //     { idx_min=node_idx; break;}
                min_dist = nodes_[node_idx].get_dist_to_config(config);
                idx_min = node_idx;
                // std::cout << "min_dist " << min_dist << ", idx_min " << idx_min << std::endl;
            }
        }
        std::cout << std::endl;
        return idx_min;
    }

    void insert_node(const Config &new_config)
    {
        if (nodes_.size()==0)
        {
            Node start_node(start_config_, -1);
            start_node.numofDOFs_ = numofDOFs_;
            nodes_.push_back(start_node);
            std::cout << "inserted first node" << std::endl;
            return;
        }

        int idx_nearest = get_nearest_index(new_config);
        std::cout << "got nearest index " << idx_nearest << std::endl; 

        // if new_config is not in collision
        if (IsValidArmConfiguration(new_config, numofDOFs_, map_, map_x_size_, map_y_size_))
        {
            // std::cout << "valid. adding to tree " << std::endl;
            // move by epsilon
            for (int joint_idx=0; joint_idx<numofDOFs_;joint_idx++)
            {
                // std::cout << " nearest " << nodes_[idx_nearest].get_config()[joint_idx];
                // std::cout << " sample " << new_config[joint_idx];
                // double nearest = nodes_[idx_nearest].get_config()[joint_idx];
                // double sample = new_config[joint_idx];
                new_config[joint_idx] = nodes_[idx_nearest].get_config()[joint_idx] +
                            (epsilon_*( new_config[joint_idx] - nodes_[idx_nearest].get_config()[joint_idx] ));// TODO weird bug with non hardcoded epsilon
                // std::cout << " to be added " << new_config[joint_idx] << std::endl;;
            }

            Node new_node(new_config, idx_nearest);
            nodes_.push_back(new_node);

            if (new_node.is_in_goal_region(goal_config_))
            {
                is_goal_reached_ = true;
                construct_path();
            }
            return;
        }

        else
        {
            // std::cout << "invalid config " << std::endl;
            return;
        }

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
        std::cout << "num_samples_ : " << num_samples_ << std::endl; 
        std::cout << "num_nodes : " << nodes_.size() << std::endl; 
        num_samples_++;

        if ((double) rand()/(double) RAND_MAX < sample_goal_bias_)
        {
            // std::cout << "sampling goal config" << std::endl;
            return goal_config_;
        }
        else
        {
            Config sample_config = new double[numofDOFs_];
            // std::cout << "sampling random" << std::endl;

            for(int i = 0; i < numofDOFs_; ++i) 
            {
                // std::cout << "i " << i << std::endl;
                sample_config[i] = (double)rand()*(2*M_PI)/(double) RAND_MAX;
                // std::cout << "sample_config[i]" << sample_config[i] << std::endl;
            }

            // std::cout << "sampling done" << std::endl;
            return sample_config;
        }

    } 

    std::vector<Config> get_path()
    {
        return path_;
    }

    void set_map(double* map, int x_size, int y_size)
    {
        map_ = map;
        map_x_size_ = x_size;
        map_y_size_ = y_size;
    }

    bool solve()
    {
        for (int iter=0; iter<num_max_iters_; iter++)
        {
            if (is_goal_reached_)
            {
                std::cout << " FOUND SOLUTION " << std::endl;
                std::cout << " FOUND SOLUTION " << std::endl;
                std::cout << " FOUND SOLUTION " << std::endl;
                std::cout << " FOUND SOLUTION " << std::endl;
                std::cout << " FOUND SOLUTION " << std::endl;
                std::cout << " FOUND SOLUTION " << std::endl;
                std::cout << " FOUND SOLUTION " << std::endl;
                return true;
            }

            Config sample_config = get_random_sample();
            insert_node(sample_config);
        }

        return false;
    }

    Config start_config_;
    Config goal_config_;
    std::vector<Node> nodes_;
    std::vector<Config> path_;
    int num_samples_;
    int num_max_iters_;
    bool is_goal_reached_;
    int numofDOFs_;
    double* map_;
    int map_x_size_;
    int map_y_size_;
    double sample_goal_bias_;
    double epsilon_;
};

static void planner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    //no plan by default
    *plan = NULL;
    *planlength = 0;
    double epsilon = 0.1;
    int num_max_iters = (int)1e2;
    double sample_goal_bias = 0.3;
    // static
    Node temp(armstart_anglesV_rad);
    temp.set_numofDOFs(numofDOFs);

    // // Config start_config(armstart_anglesV_rad, -1);
    // // Config goal_config(armgoal_anglesV_rad, -1);
    Tree rrt(armstart_anglesV_rad, armgoal_anglesV_rad, num_max_iters, epsilon, numofDOFs, sample_goal_bias);
    rrt.set_map(map, x_size, y_size);
    bool got_path = rrt.solve();
    if (got_path)
    {
        std::vector<Config> path_vector = rrt.get_path();
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
    //if (planner_id == RRT)
    //{
    //    plannerRRT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    //}
    
    //dummy planner which only computes interpolated path
    planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    
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





