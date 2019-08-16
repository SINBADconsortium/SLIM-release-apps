function [f,g,h] = LS_misfit(x,y,i_src,i_freq,C)
% LS_MISFIT - Standard least squares misfit
%
% Curt Da Silva, 2016
% 
% Usage:
%  [f,g,h] = LS_misfit(x,y,i_src,i_freq,{C});
% 
% Input:
%   x      - predicted data
%   y      - actual data 
%   i_src  - source index (only used if C is specified,optional) 
%   i_freq - frequency index (only used if C is specified,optional)
%   C      - offset mask, nrec x nsrc matrix (default: ones(nrec,nsrc))
%
% Output:
%   f      - misfit function
%   g      - misfit gradient
%   h      - misfit hessian (spot operator)
%
    if exist('C','var')==0, C = ones(size(x));
    else C = C(:,i_src); end
    r = C.*(x-y);
    f = 0.5*norm(r,'fro')^2;
    g = r;
    h = opDiag_swp(vec(C));
end
