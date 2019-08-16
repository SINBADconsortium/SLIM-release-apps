## There are four .m files:
  
1. run_TimeJitAcq_Deblending.m 
   - main script that simulates the 2-D time-jittered blended marine acquisition scenario, and executes the deblending algorithm.

2. TimeJitAcq_params.m
   - script to set the parameters used in script 1.

3. setpath.m
   - script to set the required paths (used inside script 2).

4. view_results.m
   - script to view the experimental results (output from script 1).

NOTE: the function that executes the deblending algorithm, TimeJitAcq_Deblend.m, can be found under the /misc_funcs directory.


## Running
   
Run the scripts in the following order:

- first, set the parameters in script 3
  [NOTE: you do NOT need to 'run' this script. Just use it to set the parameters. This script is called inside the main script (1).]
- second, run the main script (1)
- third, view the results using script 5


## System requirements and runtime

This example requires about 1.0 GB RAM. For parallel computation, use 10 workers (for e.g., configuration: 1 node x 10 workers).
Runtime: approximately 12 hours.


## Data adaptation

To run the algorithm on your own data, edit the parameters in the TimeJitAcq_params.m accordingly.

