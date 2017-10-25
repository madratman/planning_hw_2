RRT, RRT-Connect, RRT* are working. 
PRM is implemented but has a bug unfortunately which I am not able to figure out. 

"results_benchmarking/" has the raw results for each planner. 

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