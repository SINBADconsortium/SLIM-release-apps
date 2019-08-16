%% Fast least-squares imaging with source estimation using multiples
% This applications is available only in the software release for
% members of SINBAD consortium.
%
% This software provides an algorithm to do fast least-squares
% migration while estimating the source wavelet on the fly by variable
% projection. Using multiples further improves image quality and
% enables the retrieval of the true-amplitude image and source
% wavelet. The inversion process is sped up by using source/frequency
% subsampling and rerandomization.
%
% Author: Ning Tu (tning@eos.ubc.ca)
% Date: Dec/17/2015

%% Downloading & Dependencies
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM software release> under
% |/applications/Imaging/L1MIGRATIONwVP|.
%
% The code has been tested with _Matlab R2013b_ and requires the Parallel
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
% * _tools/solvers/GenSPGL1_ - Generalized SPGL1
%
% If you want to use your own modules to do modelling, please contact
% the author.

%% Running & Parallelism
% All the examples and results are produced by the scripts found in
% this software release under |/applications/Imaging/L1MIGRATIONwVP/examples/|. 
% Start matlab from |/applications/Imaging/L1MIGRATIONwVP| to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file
% enclosed inside the folder for each set of examples.
% 
% The scripts can be run in serial mode while parallel mode is advised.
% 
% Read more about <https://slim.gatech.edu/SoftwareDemos/tools/ SLIM's approach to parallel computing in Matlab>.

%% Examples and results
% <example.html Examples and results are shown here.>
% Scripts to reproduce the results can be found under
% |/applications/Imaging/L1MIGRATIONwVP/examples|.

%% References
% <https://slim.gatech.edu/content/source-estimation-surface-related-multiples–-fast-ambiguity-resolved-seismic-imaging
% [1]>. Ning Tu, Aleksandr Aravkin, Tristan van Leeuwen, Tim Lin, and Felix J. 
% Herrmann, "Source estimation with surface-related multiples—fast 
% ambiguity-resolved seismic imaging", submitted to Geophysical Journal 
% International. 2015
%
% <https://slim.gatech.edu/content/fast-least-squares-migration-multiples-and-source-estimation
% [2]>. Ning Tu, Aleksandr Y. Aravkin, Tristan van Leeuwen, and Felix
% J. Herrmann, “Fast least-squares migration with multiples and source
% estimation”, EAGE. 2013.
%
% <https://slim.gatech.edu/content/sparse-seismic-imaging-using-variable-projection
% [3]>. Aleksandr Y. Aravkin, Tristan van Leeuwen, and Ning Tu,
% “Sparse seismic imaging using variable projection”, ICASSP. 2012
%


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support. The
% author is also receives a scholarship from China Scholarship
% Council. The synthetic Compass model is by courtesy of BG Group.
