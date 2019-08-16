## There are four .m files:
  
1. SourceSep_HSS.m 
   - main script that simulates the 2-D over/under blended marine acquisition scenario, and executes the source separation algorithm.

2. SourceSep_params.m
   - script to set the parameters used in 1.

3. setpath.m
   - script to set the required paths (used inside 2).

4. view_results.m
   - script to view the experimental results (output from 1).


## Running
   
Run the scripts in the following order:

- first, set the parameters in script 2
  [NOTE: you do NOT need to 'run' this script. Just use it to set the parameters. This script is called inside the main script (1).]
- second, run the main script (1)
- third, view the results using script 4


## System requirements and runtime

This example requires about 2 GB RAM. Runtime: 24 hours (approx.)


## Data adaptation

To run the algorithm on your own data, edit the parameters in the
SourceSep_params.m and SourceSep_HSS.m scripts accordingly.

