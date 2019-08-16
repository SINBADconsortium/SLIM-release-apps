function [V] = InterpModel(m,n)
%% function V = InterpModel(M,n)
% 2D interpolate a velocity model M to size n
nt = size(m);
xt = [0 : 1 : nt(2)-1];
zt = [0 : 1 : nt(1)-1];

x  = linspace(xt(1),xt(end),n(2));
z  = linspace(zt(1),zt(end),n(1));

Sx = opLInterp1D(xt(:),x(:));
Sz = opLInterp1D(zt(:),z(:));
S   = opKron(Sx,Sz);
V   = reshape(S*m(:),n);
