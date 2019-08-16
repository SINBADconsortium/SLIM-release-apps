function A = LinInterp1D(xin,xout)
% 1D Linear interpolation matrix at arbitrary locations
% 
% Author: Curt Da Silva, 2015
%
% Usage:
%   A = LinInterp1D(xin,xout);
% 
% Input:
%   xin - input grid
%   xout - output grid (must be contained in [min(xin),max(xin)])
%
% Output:
%   A   - length(xout) x length(xin) sparse matrix
    
    xin = sort(vec(xin),'ascend');
    xout = sort(vec(xout),'ascend');
    nin = length(xin); 
    nout = length(xout);
    I = nearestnbr1d(xin,xout); 
    I = vec(I);
    I(xin(I)>xout) = I(xin(I)>xout)-1;
    xin = [xin;1e12];
    
    b = (xout-xin(I))./(xin(I+1)-xin(I));
    a = (xin(I+1)-xout)./(xin(I+1)-xin(I));
    A = sparse(vec(1:nout),I,vec(a),nout,nin) + sparse(vec(1:nout),min(I+1,nin),vec(b),nout,nin);
    
end