function [rho, x] = CalDist(smp, k)

if nargin < 2
        k =100;
 end
 
 for i = 1: size(smp,1)
        [N X]    = hist(smp(i,:),k);
%        rho(i,:) = N / max(N);
        rho(i,:) = N / size(smp,2);
        x(i,:)      = X; 
        rho(i,:) = rho(i,:)/(x(i,2)-x(i,1));
 end

