function x = ifnt1d(y)
% Fast 1-D Inverse Noiselet Transform
%
% Description :
%   This algorithm is implemented in N log2 N additions and subtractions.
% In: 
%   y : Initial data transformed by fnt1d. Length should be an integer power of 2.
% Out:
%   x : the noiselet inverse tranform of y
% Example:
%   x = rand(32,1);
%   y = fnt1d(x);
%   nx = ifnt1d(y);
%   max(abs(x-nx)) % about 1e-16
% See also:
%   ifnt1d
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

x = bitrevorder(y);

if isvector(y)
    N = length(y);
    x = reshape(x, [N 1]);
else
    N = size(y, 1);
end

ndigit = floor(log2(N));

cp = 1 + i; 
cm = 1 - i;

ind = uint32(0:(N-1));
pos = (mod(ind,2) == 1);

for i1 = 1:ndigit   
    pos = [pos(1:2:end) pos(2:2:end)];

    ni = ind(~pos) + 1;
    nj = ind(pos) + 1;
        
    temp1 = x(ni, :); temp2 = x(nj, :); 
    x(ni, :) = cp*temp1 + cm*temp2;
    x(nj, :) = cm*temp1 + cp*temp2;
end

x = inv(N)*x(end:-1:1, :); %Delete this line for inverse transform
