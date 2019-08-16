function y = fnt1d(x)
% Fast 1-D Noiselet Transform
%
% Description :
%   This algorithm is implemented in N log2 N additions and subtractions ().
% In: 
%   x : Data sequence. Length should be an integer power of 2. If x is an
%   array, fnt1d is performed on the columns (with the same behaviour
%   than fft on  arrays)
% Out:
%   y : the noiselet tranform of x
% Example:
%   x = rand(32,1);
%   y = fnt1d(x);
% See also:
%   ifnt1d, fnt1d, fnt2d, ifnt2d
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

y = bitrevorder(x);

if isvector(x)
    N = length(x);
    y = reshape(y, [N 1]);
else
    N = size(x, 1);
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
        
    temp1 = y(ni,:); temp2 = y(nj,:); 
    y(ni,:) = cp*temp1 + cm*temp2;
    y(nj,:) = cm*temp1 + cp*temp2;
end

y = inv(N)*y; 


