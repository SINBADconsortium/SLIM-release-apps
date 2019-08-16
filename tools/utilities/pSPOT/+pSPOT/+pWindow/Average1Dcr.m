function [ t ] = Average1Dcr( i, n, h )
%Average1Dcr 1D taper using .5 - for outer windows
%
%   [ T ] = Average1Dcr( I, N, H )
%
%   INPUT:
%      I = element index (1:N)
%      N = window length
%      H = half of the overlap's size
%   OUTPUT:
%      T = taper's value at I
%

    if i<=n-2*h
        x=n;
    else
        x=n-i+1;
    end
    b=2*h;
    if x<n;
        t=.5;
    else
        t=1;
    end
%    fprintf('%2d %2d %2d %f\n',b,i,x,t);
end

