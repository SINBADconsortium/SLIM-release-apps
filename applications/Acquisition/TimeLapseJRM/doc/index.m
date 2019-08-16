%% Joint recovery method for time-lapse seismic data
%
% This application is available only in the software release for members of SINBAD consortium.
% 
% This package is an application for time-lapse data acquired with randomized sampling schemes, and a joint recovery method based on sparse inversion (via L1 minimization). In particular, we show the performance of the method for: (i) data with random missing shots, and (ii) data simulated for a time-jittered, blended marine acquisition.
% 
%
% Authors: Felix Oghenekohwo (foghenekohwo@eos.ubc.ca), Haneet Wason (hwason@eos.ubc.ca)
%
% Date: June, 2014


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/Acquisition/TimeLapseJRM|.
%
% The code has been tested with _Matlab R2013a_.
%
% This code uses the following packages, also found in the |tools| part of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _functions/misc_ - miscellaneous functions.
% * _solvers/SPGL1-SLIM_ - SLIM version of SPGL1 solver.
% * _transforms/CurveLab-2.1.2-SLIM_ - curvelet transform functions.


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in this software release under |applications/Acquisition/TimeLapseJRM/|. Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed with the code. The scripts can be run in serial mode.
%
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.


%% Examples and results
%
% <results.html An example of recovery from data with missing shots is shown here.>
% 
% <example_marine.html An example of recovery from time-jittered, blended data is shown here.> 


%% References
%
% <https://slim.gatech.edu/content/time-lapse-seismic-without-repetition-reaping-benefits-randomized-sampling-and-joint-recover [1]> Felix Oghenekohwo, Ernie Esser, and Felix J. Herrmann [2014]. Time-lapse seismic without repetition: reaping the benefits from randomized sampling and joint recovery. Presented at the 76th EAGE Conference and Exhibition.
%
% <https://slim.gatech.edu/content/randomization-and-repeatability-time-lapse-marine-acquisition-0 [2]> Haneet Wason, Felix Oghenekohwo, and Felix J. Herrmann [2014]. Randomization and repeatability in time-lapse marine acquisition. Presented at the EAGE Workshop on Land and Ocean Bottom; Broadband Full Azimuth Seismic Surveys, Spain.
%
% <https://slim.gatech.edu/content/randomization-and-repeatability-time-lapse-marine-acquisition [3]> Haneet Wason, Felix Oghenekohwo, and Felix J. Herrmann [2014]. Randomization and repeatability in time-lapse marine acquisition. To be presented at the SEG Conference.


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.

