function X = ifnt2d(Y)
% Fast 2-D Inverse Noiselet Transform
%
% Description :
%   This algorithm is implemented in N^2 log2^2 N additions and subtractions.
% In: 
%   Y : The Noiselet Transform of an Image produced by fnt2d() 
% Out:
%   Y : the 2-D noiselet tranform of it
% Example:
%   X = rand(32,32);
%   nX = ifnt2d(fnt2d(X));
%   max(abs(X(:)-nX(:))) % About 1e-16
% See also:
%   fnt2d, fnt1d, ifnt1d
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


X = ifnt1d(ifnt1d(Y).').';