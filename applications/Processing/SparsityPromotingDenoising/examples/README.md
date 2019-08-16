## There are four .m files:
  
1. getCurveletCoeff.m 
   - script that computes the curvelet coefficients of the input seismic data.

2. denoiseSeismicData.m
   - script that runs the denoising algorithm.

3. setDenoiseParams.m
   - script to set the parameters used in scripts 1 and 2.

4. setpath.m
   - script to set the required paths (used inside script 3).


## Running
   
Run the scripts in the following order:

- first, set the parameters in script 3.
  [NOTE: you do NOT need to 'run' this script. Just use it to set the parameters. This script is called inside scripts 1 and 2.]
- second, run script 1 to get the curvelet coefficients, which will be used for denoising in script 2.
- third, run script 2 to denoise data.


## System requirements and runtime

This example runs in serial mode. \
**NOTE** - the parameters (in script 3, setDenoiseParams.m) for this example are set for dynamic geometry with 'WRAP' curvelets.

### Memory requirements 

Most of the memory required is when running script 1 (getCurveletCoeff.m). Approximate memory required is about: \
1. Dynamic geometry (DG): 'WRAP' (wrapping) curvelets - 6 GB; 'ME' (mirror extended) curvelets - 7 GB \
2. Static geometry (SG): 'WRAP' curvelets - 15.5 GB; 'ME' curvelets - 22 GB  

For script 2 (denoiseSeismicData.m), the approximate memory required is about: \
1. 'WRAP' curvelets: DG and SG - 2 GB \
2. 'ME' curvelets: DG and SG - 5 GB

### Runtime

Most of the computation time is spent in script 1 (getCurveletCoeff.m). \
1. Dynamic geometry (DG): 'WRAP' curvelets - 1 hour; 'ME' curvelets - 2.4 hours  
2. Static geometry (SG): 'WRAP' curvelets - 6.1 hours; 'ME' curvelets - 26 hours  

Script 2 is run in interactive mode as the user tries different thresholding values to denoise input data and view results (interactively) until the desired output is achieved.


## Data adaptation

To run the algorithm on your own data, edit the parameters in the setDenoiseParams.m and denoiseSeismicData.m accordingly. \
**NOTE** - this example denoises frequency slices, however, the scripts can be adapted to denoise time slices.

