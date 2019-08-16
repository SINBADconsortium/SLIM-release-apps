%% Time-domain Gauss-Newton full-waveform inversion (applied to Chevron 2014 benchmark dataset)
%
% This application is available only in the software release for members of SINBAD consortium.
%
% This software release includes an parallel framework in Matlab for time-domain 
% Gauss-Newton full-waveform inversion [1]. This application package is an extension 
% of our previous work which is conducted in the frequency domain [2,3]. In this package, we
% only carried out Gauss-Newton FWI with randomized dimensionality reduction techniques.
% Our future work will be time domain FWI with modified Gauss-Newton techniques [2,3].
%
% Author: Dr. Xiang Li
%
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date  : July, 2015

%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware/ SLIM software release> under
% |/applications/WaveformInversion/2DBasicTimeStepping|.
%
% The code has been tested with _Matlab R2013a_ and requires the Parallel
% Computing Toolbox.
%
% To run this code, please download Chevron 2014 SEG workshop benchmark dataset from 
% <https://s3.amazonaws.com/open.source.geoscience/open_data/seg_workshop_fwi_2014/seg_workshop_fwi_2014.html>
% You can also find DATA LICENSE AGREEMENT information with the above link. Please
% use this dataset according to thire requirement presented in the above link.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware/ SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/TimeModeling - 2D/3D time-stepping modeling kernel
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.


%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |applications/WaveformInversion/2DBasicTimeStepping/|. 
% Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instrictions in the README file enclosed
% with the code.
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
%
% Scripts to generate results of examples can be found in the scripts directory. <examples_TSGNFWI.html
% The results are shown here>.
%

%% References
%
% <https://slim.gatech.edu/content/2d3d-time-stepping-imaging-and-inversion> [1]. Xiang Li - 2D/3D time-stepping for imaging and inversion, SINBAD Fall consortium talks. 2014
%
% <https://slim.gatech.edu/node/6390 [2]> Felix J. Herrmann, Xiang Li, Aleksandr Y. Aravkin, and Tristan van Leeuwen, A modified, sparsity promoting, Gauss-Newton algorithm for seismic waveform inversion, in Proc. SPIE, 2011, vol. 2011.
%
% <https://slim.gatech.edu/node/6621 [3]> Xiang Li, Aleksandr Y. Aravkin, Tristan van Leeuwen, and Felix J. Herrmann, Fast randomized full-waveform inversion with compressive sensing. 2011. Geophysics, accepted.
%
