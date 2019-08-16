## There are two .m files:
  
1. JRM_OneRecvGather.m 
   - main script that simulates the data acquisition with missing shots, 
   - seen on one common receiver gather
   - showing different realizations for baseline/monitor data
   - and executes the joint recovery method
   - to reconstruct the fully sampled timelapse wavefields

2. setpath.m
   - script to set the required paths (used inside 1).

3. JRM_mkfigs.m
   - script to view the experimental results (output from 1).


## Running
   
Run the scripts in the following order:
- first (optional), change the subsampling parameters in script 1 (line 31) 
- second, run the main script (1)
- third, view the results using script 3


## Data adaptation

To run the algorithm on your own data, implementation procedure is 
in progress.

