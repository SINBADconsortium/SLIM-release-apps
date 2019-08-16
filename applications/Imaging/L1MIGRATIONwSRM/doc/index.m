%% Fast imaging with surface-related multiples
% This applications is available only in the software release for
% members of SINBAD consortium.
%
% This software provides an algorithm to image from the entire upgoing
% wavefield by sparse inversion. The inversion process is sped up by
% using source/frequency subsampling and renewal.
%
% Author: Ning Tu (tning@eos.ubc.ca)
% Date: Dec/17/2015

%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |/applications/Imaging/L1MIGRATIONwSRM|.
%
% The code has been tested with _Matlab R2012b_ and requires the Parallel
% Computing Toolbox.
%
% This code uses the following packages, also found in the |tools| part
% of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/pSPOT_ - parallel extension of SPOT.
% * _algorithms/2DFreqModeling_ - 2D constant-density acoustic modeling
% * _tools/algorithms/REPSI_ - Robust EPSI algorithm
% * _operators/misc_ - Misc. operators such as interpolation, smoothing and splines
% * _functions/misc_ - Misc. functions.
%
% If you want to use your own modules to do modelling or multiple
% prediction, please contact the author.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/Imaging/L1MIGRATIONwSRM/examples/|. 
% Start matlab from |/applications/Imaging/L1MIGRATIONwSRM| to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file
% enclosed inside the folder for each set of examples.
% 
% The scripts can be run in serial mode but parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
% <example.html Examples and results are shown here.>
% Scripts to reproduce the results can be found under
% |/applications/Imaging/L1MIGRATIONwSRM/examples|.

%% References
% <https://slim.gatech.edu/content/fast-imaging-surface-related-multiples-sparse-inversion [1]>. Tu, Ning and Felix J. Herrmann,
% [2012] Fast imaging with surface-related multiples by sparse inversion.
% Geophysical Journal International, vol 201, 304-317. See the "References"
% section of this paper for a detailed list of references.

%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support. The author is
% also financially supported by China Scholarship Council. The salt dome model
% is by courtesy of Eric Verschuur. The Sigsbee 2B model is by courtesy of the
% SMAART JV.
