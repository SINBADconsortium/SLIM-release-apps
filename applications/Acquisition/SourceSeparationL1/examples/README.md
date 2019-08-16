## There are four .m files:
  
1. SourceSepL1.m 
   - main script that executes the source separation algorithm.

2. SourceSepL1_params.m
   - script to set the parameters used in 1.

3. setpath.m
   - script to set the required paths (used inside 2).

4. convert_results_freqtotime.m
   - script to convert the experimental results (output from 1) from frequency domain to time domain.


## Running
   
Run the scripts in the following order:

- first, set the parameters in script 2
  [NOTE: you do NOT need to 'run' this script. Just use it to set the parameters. This script is called inside the main script (1).]
- second, run the main script (1)
- once the experiment has finished, return to this directory to view the results using script 4


## System requirements and runtime

See the script 1 (SourceSepL1.m) to know how the code can be edited to run in parallel (recommended).
This example requires about 7 GB RAM. Runtime: serial - 45 hours (approx.); parallel (10 sets) - 4.5 hours (approx.).


## Data adaptation

To run the algorithm on your own data, edit the parameters in the
SourceSepL1_params.m and SourceSepL1.m scripts accordingly.

