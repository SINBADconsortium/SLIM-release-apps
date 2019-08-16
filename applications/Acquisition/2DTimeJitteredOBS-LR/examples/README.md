## There are four .m files:
  
1. jitacq_deblending_LR.m 
     - main script that simulates the 2D time-jittered marine acquisition scenario, and 
       executes the deblending algorithm.

2. jitacq_deblending_LR_params.m 
     - script to set the time-jittered acquisition parameters used in 1.

3. setpath.m
     - script to set the required paths (used inside 2)

4. view_results.m
     - script to view the experimental results (output from 1)


## RUNNING

Run the scripts in the following order:

- first, set the time-jittered acquisition parameters (script 2)
  [NOTE: you do NOT need to 'run' this script. Just use it to set the parameters.
         This script is called inside the main script (1).]
- second, run the main script (1) 
- third, view the results using script 4


## SYSTEM REQUIREMENTS AND RUNTIME

This example requires about 16 GB RAM. This version runs in serial mode.
     
Runtime: 2 hours (approx.)


## Testing on your own data
   
   To run the algorithm on your own data, run the scripts in the following order:
     - first, set the time-jittered acquisition parameters (script 2)
       [NOTE: you do NOT need to 'run' this script. Just use it to set the parameters.
              This script is called inside the main script (1).]
     - second, run the main script (1) 
     - third, view the results using script 3

