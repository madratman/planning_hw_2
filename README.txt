RRT, RRT-Connect, RRT* are working. 
PRM is implemented but has a bug unfortunately which I am not able to figure out. 

The implementations are pretty much from the snapshots of the algorithms from the slides. 

For RRT*, for the radius calculation:
 - gamma used was 100
 - epsilon (epsilon_rrt_star_ below) was PI/4 
 (from piazza discussion and paper, the slides are wrong as it's not number of vertices.
 - delta_ is volume of R5 hyperball, aka pow(M_PI, 2.5) / gamma_function(3.5)
radius_rrt_star_ = std::min(pow( gamma_ / delta_ * log(nodes_.size()) / nodes_.size() , 1.0/(double)numofDOFs_), epsilon_rrt_star_);

For collision checking, I interpolate 20 times between current config and sample config, and check each intermediate config for collision. 
For the extend function in RRT-Connect, I again extend by interpolation step size of 20.  
In PRMs, I keep on adding nodes in a while(true) loop till the graph is fully connected.  

I mostly maintain a vector of a Node class, and keep track of indices. Should have used pointers in hindsight, but I guess effectively it's the same thing, but readable and easier to write code. 

Plan quality is just cumulative sum of (minimum) joint angle differences over the path returned.

random_samples.txt has 40 valid joint configs. I simply loop over them in (start_config = (2i-1)th, goal_config(2i)th fashion, as can be seen in rrt.m. 

The code can be run as asked in the PDF:
runtest('map2.txt',startQ, goalQ, 0/1/2);

Also, this is definitelythe most useful homework assignments I have ever done! Thanks!

The "results_benchmarking/" folder has the raw results for each planner. 

+-------------+--------------------+--------------------+-------------------+
|   Planner   | Mean No of Samples | Mean planning time | Mean path quality |
+-------------+--------------------+--------------------+-------------------+
|     RRT     |     1107.750000    |      0.027118      |      5.748375     |
+-------------+--------------------+--------------------+-------------------+
| RRT-Connect |     680.700000     |      0.007917      |      5.603374     |
+-------------+--------------------+--------------------+-------------------+
|   RRT-Star  |     3202.700000    |      0.319238      |      6.384765     |
+-------------+--------------------+--------------------+-------------------+
Table generator credits : http://www.tablesgenerator.com/text_tables