### Implementation of Sampling Based Motion Planning Algorithms for a 5 DoF planar arm In C++

Implementation of Sampling Based Motion Planning Algorithms as HW2 solution of for course [16-782 Planning and Decision-making in Robotics](http://www.cs.cmu.edu/~maxim/classes/robotplanning_grad/) at CMU:
- RRT
- RRT*
- RRT-Connect
- PRM 
for a planar, 5 DoF arm 

HW assignment available [here](http://www.cs.cmu.edu/~maxim/classes/robotplanning_grad/hw/hw2_16782_fall17.pdf)   
Will be making this agnostic to state space soon!

### How to run

### Results.
Each planner was run with random start and goal configurations (found in `random_samples.txt` (odd lines were chosen as start, even as end config of joint angles).

Path quality is euclidean difference between each "waypoint" in joint space

|   Planner   | Mean No of Samples | Mean planning time | Mean path quality |
|:-----------:|:------------------:|:------------------:|:-----------------:|
|     RRT     |     1107.750000    |      0.027118      |      5.748375     |
| RRT-Connect |     680.700000     |      0.007917      |      5.603374     |
|   RRT-Star  |     3202.700000    |      0.319238      |      6.384765     |
