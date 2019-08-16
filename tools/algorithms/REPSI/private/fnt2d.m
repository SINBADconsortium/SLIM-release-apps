function Y = fnt2d(X)
% Fast 2-D Noiselet Transform
%
% Description :
%   This algorithm is implemented in N^2 log2^2 N additions and subtractions.
% In: 
%   X : The image sequence. Image dimensions should be an integer power of 2. 
% Out:
%   Y : the 2-D noiselet tranform of it
% Example:
%   X = rand(32,32);
%   Y = fnt2d(X);
%   max(ifnt1d(fnt1d(x)))
% See also:
%   ifnt2d, fnt1d, ifnt12
% Author:
%   L. Jacques (LTS2/EPFL, 2008)
%
% References :
% * This code is just a slight variation of the Fast (Cooley-Tukey type)
%   (1D and 2D) Hadamard Transform, by Gylson Thomas :
%     fwht1d.zip: 
%          http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6879 (fwht1d.zip)
%     fwht2d.zip:
%          http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6882
% * R. Coifman, F. Geshwind, and Y. Meyer. Noiselets. Appl. Comput. Harmon. Anal., 10(1):27--44, 2001.

Y = fnt1d(fnt1d(X).').';