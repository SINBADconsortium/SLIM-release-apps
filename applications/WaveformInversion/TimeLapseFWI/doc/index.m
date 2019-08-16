%% Full-Waveform Inversion for time-lapse seismic data
%
% This application is available only in the software release for members of SINBAD consortium.
% 
% This package is an application for full-waveform inversion of time-lapse seismic data sets using a modified Gauss-Newton inversion framework. Ideas from distributed compressive sensing and stochastic optimization are exploited to create improved time-lapse models and speed-up the inversion, respectively.
% 
%
% Author: Felix Oghenekohwo (foghenekohwo@eos.ubc.ca)
%
% Date: February, 2015


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/WaveformInversion/TimeLapseFWI|.
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
% All the examples and results are produced by the scripts found in this software release under |applications/WaveformInversion/TimeLapseFWI/|. Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed with the code. The scripts can be run in serial mode.
%
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.


%% Examples and results
%
% <example_timelapse_BG.html An example of our joint FWI on the BG time-lapse model is shown here.> 


%% References
%
% <https://slim.gatech.edu/content/time-lapse-seismic-without-repetition-reaping-benefits-randomized-sampling-and-joint-recover [1]> Felix Oghenekohwo, Ernie Esser, and Felix J. Herrmann [2014]. Time-lapse seismic without repetition: reaping the benefits from randomized sampling and joint recovery. Presented at the 76th EAGE Conference and Exhibition.
%
% <https://slim.gatech.edu/node/6390 [2]> Felix J. Herrmann, Xiang Li, Aleksandr Y. Aravkin, and Tristan van Leeuwen, A modified, sparsity promoting, Gauss-Newton algorithm for seismic waveform inversion, in Proc. SPIE, 2011, vol. 2011.
%
% <https://slim.gatech.edu/node/6621 [3]> Xiang Li, Aleksandr Y. Aravkin, Tristan van Leeuwen, and Felix J. Herrmann, Fast randomized full-waveform inversion with compressive sensing. 2011. Geophysics.



%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support. The synthetic time-lapse model was provided by BG Group.

