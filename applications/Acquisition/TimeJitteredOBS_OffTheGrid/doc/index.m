%% Time-jittered blended marine acquisition on non-uniform spatial grids
%
% This application is available only in the software release for members of SINBAD consortium.
% 
% This package contains a MATLAB implementation of a 2-D time-jittered blended marine acquisition scheme on non-uniform spatial (source) grid, and a deblending algorithm based on sparse inversion via one-norm-minimization incorporating the non-equispaced fast discrete curvelet transform.
% 
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: August, 2015


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |applications/Acquisition/TimeJitteredOBS_OffTheGrid|.
%
% The code has been tested with _Matlab R2014b_.
%
% This code uses the following packages, also found in the |tools| part of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _functions/misc_ - miscellaneous functions.
% * _solvers/SPGL1-SLIM_ - SLIM version of SPGL1 solver.
% * _transforms/CurveLab-2.1.2-SLIM_ - curvelet transform functions.


%% Functions
%
% Some functions specific to this package can be found in the |applications/Acquisition/TimeJitteredOBS_OffTheGrid/misc_funcs| directory.


%% Running & Parallelism
%
% All the examples and results are produced by the scripts found in this software release under |applications/Acquisition/TimeJitteredOBS_OffTheGrid|. Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed with the code.
%
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.


%% Examples and results
%
% <example.html Examples and results are shown here.> Scripts to reproduce the results can be found under |applications/Acquisition/TimeJitteredOBS_OffTheGrid/examples|. 


%% References
%
% <https://slim.gatech.edu/content/time-jittered-ocean-bottom-seismic-acquisition [1]> Haneet Wason and Felix J. Herrmann [2013]. Time-jittered ocean bottom seismic acquisition. SEG Technical Program Expanded Abstracts, vol. 32, pp. 1-6.
%
% <https://slim.gatech.edu/content/nonequispaced-curvelet-transform-seismic-data-reconstruction-sparsity-promoting-approach [2]> Gilles Hennenfent, Lloyd Fenelon and Felix J. Herrmann [2010]. Nonequispaced curvelet transform for seismic data reconstruction: a sparsity-promoting approach, Geophysics, vol. 75, pp. WB203-WB210.
%
% <https://slim.gatech.edu/content/compressed-sensing-4-d-marine–-recovery-dense-time-lapse-data-subsampled-data-without-repeti [3]> Haneet Wason, Felix Oghenekohwo and Felix J. Herrmann [2015]. Presented at the 77th EAGE Conference & Exhibition (2015), doi: 10.3997/2214-4609.201413088.


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.

