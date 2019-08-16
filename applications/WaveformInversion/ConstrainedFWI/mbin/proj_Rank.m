function [out]=proj_Rank(in,r)

[U,S,V]=svd(in);
out=U(:,1:r)*S(1:r,1:r)*V(:,1:r)';