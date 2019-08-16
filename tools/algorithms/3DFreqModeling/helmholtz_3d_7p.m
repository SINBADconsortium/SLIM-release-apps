% 3D Helmholtz operator, discretized using the standard 7-point stencil with PML
%
% USE:
%   H = helmholtz_3d_7p(model,h,np,freqhz,unit)
%
% INPUT:
%   freqhz  - frequency in Hertz
%   model   - velocity model in meters per second or slowness squared; this is a 
%             *3D object* with dimension (nx,ny,nz), NOT A VECTOR!
%   h       - 3D vector containing the grid spacing in each direction
%   np      - structure containing the number of np points in each direction 
%             np.x, np.y, np.bottom and np.top (the later two stand for the z 
%             direction). For free-surface, set np.top to zero
%   unit    - String set either to "m/s" (meters per second) or "s2/m2" 
%             (slowness squared)
%
% OUTPUT:
%   H       - sparse matrix
%
% AUTHOR: Curt Da Silva, 2016
%         
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%-------------------------------------------------------------------------------
function [H,dH,ddH] = helmholtz_3d_7p(v,h,np,freq,unit)    
% size(v) <- includes pml 
   [nx,  ny,  nz] = size(v);
   N = nx*ny*nz;
   
   C = 10;
   
   pmlx = pml_func1d(nx,np(:,1),freq,C); 
   pmly = pml_func1d(ny,np(:,2),freq,C); 
   pmlz = pml_func1d(nz,np(:,3),freq,C); 
             
   [fm,dfm,ddfm] = param2wavenum(v,freq,unit);
   
   ix = vec(1:nx); iy = vec(1:ny); iz = vec(1:nz);
   
   D = @(pxM,pxN,pxP) spdiags([circshift(pxN.*pxM,-1),-pxN.*(pxP+pxM),circshift(pxN.*pxP,1)],-1:1,length(pxM),length(pxM));
   
   pxN = pmlx(ix); pxP = pmlx(ix+0.5); pxM = pmlx(ix-0.5);
   Dx = 1/h(1)^2*D(pxM,pxN,pxP);
   pyN = pmly(iy); pyP = pmly(iy+0.5); pyM = pmly(iy-0.5);
   Dy = 1/h(2)^2*D(pyM,pyN,pyP);
   pzN = pmlz(iz); pzP = pmlz(iz+0.5); pzM = pmlz(iz-0.5);
   Dz = 1/h(3)^2*D(pzM,pzN,pzP);
   M = spdiags(vec(fm),0,N,N);
   H = -(kron(speye(nz*ny),Dx) + kron(speye(nz),kron(Dy,speye(nx))) + kron(Dz,speye(nx*ny)) + M);
   dH = -spdiags(vec(dfm),0,N,N);
   ddH = -spdiags(vec(ddfm),0,N,N);
end 

function func = pml_func1d(nx,nb,freq,C)
    dist_from_int = @(x) (x <= nb(1)) .* (nb(1)-x) + (x > nx-nb(2)) .* (x-(nx-nb(2)+1));
    sigma = @(x) C * (dist_from_int(x)/max(nb)).^2;
    func = @(x) (1-1i*sigma(x)/freq).^(-1);
end