function [D,E] = getDiscreteGrad(n1,n2,h1,h2,dpw)
%n1 is the number of grid points in the vertical (depth) dimension
%n2 is the number of grid points in the lateral dimension
%h1 is the distance between grid points in the vertical (depth) dimension
%(assumed constant)
%h2 is the distance between grid points in the lateral dimension
%(assumed constant)
% dwp are the depth weights. Set to all ones for no depth weighting


%average depth weights
if (min(dpw)<=0)
    disp('problem: depth weights pm.dpw must be positive');
    return
elseif (max(dpw)>1)
    disp('problem: max depth weight should be less than one');
    return
else
    
dpw_ave = (dpw(1:end-1)+dpw(2:end))/2;

%define difference matrix D acting on vectorized model using kron
%I1 = speye(n1); %z
I2 = speye(n2); %x
Dz = spdiags(ones(n1,1)*[-1 1],0:1,n1-1,n1)/h1;
Dx = spdiags(ones(n2,1)*[1 -1],-1:0,n2,n2-1)/h2;
D = [kron(I2,spdiags(dpw_ave,0,n1-1,n1-1)*Dz) ;...
    kron(Dx',spdiags(dpw,0,n1,n1))]; %DTD is also minus discrete Laplacian
E = (D<0); %keep track of forward differences for defining vector l1 norm for isotropic TV

end