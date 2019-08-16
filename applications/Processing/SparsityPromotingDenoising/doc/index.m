%% Sparsity-promoting denoising of seismic data
%
% This application is available only in the software release for members of SINBAD consortium.
% 
% This package contains a MATLAB implementation of sparsity-promoting denoising of seismic data in the curvelet domain using one-norm minimization. For details, please see Reference [1].
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: December, 2015


%% Downloading & Dependencies
%
% The code can be found in the <https://slim.gatech.edu/consortiumsoftware SLIM sofware release> under
% |/applications/Processing/SparsityPromotingDenoising/|.
%
% The code has been tested with _Matlab R2014b_.
%
% This code uses the following packages, also found in the |tools| part of the <https://slim.gatech.edu/consortiumsoftware SLIM software release>.
%
% * _utilities/SPOT_ - object oriented framework for matrix-free linear algebra.
% * _utilities/SegyMAT_ - Matlab/Octave toolbox for reading and writing SEG-Y formatted files.
% * _solvers/SPGL1-SLIM_ - SLIM version of SPGL1 solver.
% * _transforms/CurveLab-2.1.2-SLIM_ - curvelet transform functions.


%% Running 
%
% All the examples and results are produced by the scripts found in this software release under |applications/Processing/SparsityPromotingDenoising/|. Start matlab from that directory or run |startup| in that directory to add the appropriate paths.
%
% To run the scripts follow the instructions in the README file enclosed with the code. The scripts are run in serial mode.


%% Examples and results
%
% <example.html Examples of denoising frequency slices (extracted from a 3D seismic data cube) are shown here.> 


%% References
%
% <https://slim.gatech.edu/research/processing [1]> SLIM's research webpage on processing. See the 'Denoising' section (<https://slim.gatech.edu/research/processing#denoising-via-hybrid-support-identification-and-projection>).
%
% <https://slim.gatech.edu/content/frugal-full-waveform-inversion-theory-practical-algorithm [2]> Felix J. Herrmann, Andrew J. Calvert, Ian Hanlon, Mostafa Javanmehri, Rajiv Kumar, Tristan van Leeuwen, Xiang Li, Brendan Smithyman, Eric Takam Takougang, and Haneet Wason [2013].


%% Acknowledgements
% Thanks to our sponsors and NSERC for their financial support.

