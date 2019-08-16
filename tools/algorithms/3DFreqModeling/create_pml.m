% PML function for the staggered grid 27-points stencil.
% [Operto at al. 2007. Geophysics 72(5), SM195]
%
% USE:
%   [xix, xiy, xiz] = create_pml(nx, ny, nz, np)
%
% INPUT:
%   nx      - number of points in the x direction
%   ny      - number of points in the y direction
%   nz      - number of points in the z direction
%   np      - structure containing the number of pml points in each direction
%             np.x, np.y, np.bottom and np.top
%
% OUTPUT:
%   xix     - two vectors containing the values of the ξ function in the 
%             direction x.
%   xiy     - same as above, y direction
%   xiz     - same as above, z direction
%
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: January, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%-----------------------------------------------------------------------------
function [xix, xiy, xiz] = create_pml(nx, ny, nz, np) %#codegen
   % Define the "scaled step". If the domain would be the unitary cube, s_hx would be hx (the 
   % distance between each point in the grid, the same used for discretizing the operator later).
   s_hx=1/(nx+1);
   s_hy=1/(ny+1);
   s_hz=1/(nz+1);

   % Initialize with zeroes
   xix = zeros(nx,2);
   xiy = zeros(ny,2);
   xiz = zeros(nz,2);
         
   % Scaled PML width
   if isa(np,'struct')
               
   else
      npx_top = np(1,1); npx_bot = np(2,1); 
      npy_top = np(1,2); npy_bot = np(2,2);
      npz_top = np(1,3); npz_bot = np(2,3);
   end
   Lx_top = npx_top/nx; Lx_bot = npx_bot/nx;
   Ly_top = npy_top/ny; Ly_bot = npy_bot/ny;
   Lz_top = npz_top/nz; Lz_bot = npz_bot/nz;
   
   gamma = zeros(nx+2,1);
   xix = complex(zeros(nx,2),zeros(nx,2));
   xiy = complex(zeros(ny,2),zeros(ny,2));
   xiz = complex(zeros(nz,2),zeros(nz,2));
   %------------------------------------------------------------------------------------------------
   % Computation of ξ
   %------------------------------------------------------------------------------------------------
   % The value of ξ is computed in a weighted fashion. Let φ(u) = 1 + iγ(u). Then:
   % 
   % ξ_1(u) = ____________2___________       and       ξ_2(u) = _____________2______________
   %           φ(u+1)×( φ(u+1) + φ(u) )                         φ(u+1)×( φ(u+1) + φ(u+2) )
   % 
   % - Lago
   %------------------------------------------------------------------------------------------------
   
   % For X Axis
   % Use γ to compute PML function; there are two boundaries in the
   % x axis. We have 2 of each because the grid is staggered
   %-------------------------------------------------------------------
   gamma(1:npx_bot)           = cos((pi*(0:npx_bot-1)   *s_hx) /(2*Lx_bot));
   gamma(npx_bot+1:nx+2-npx_top) = 0;
   gamma(nx+3-npx_top:nx+2)   = cos((pi*(1-(nx+2-npx_top:nx+1)*s_hx))/(2*Lx_top));
   xix(1:nx,1)=2./((1 + 1j*gamma(2:nx+1)).*((1 + 1j*gamma(2:nx+1)) + (1 + 1j*gamma(1:nx))) );
   xix(1:nx,2)=2./((1 + 1j*gamma(2:nx+1)).*((1 + 1j*gamma(2:nx+1)) + (1 + 1j*gamma(3:nx+2))) );

   % For Y Axis
   %-------------------------------------------------------------------
   gamma(1:npy_bot)           = cos((pi*(0:npy_bot-1)   *s_hy) /(2*Ly_bot));
   gamma(npy_bot+1:ny+2-npy_top) = 0;
   gamma(ny+3-npy_top:ny+2)   = cos((pi*(1-(ny+2-npy_top:ny+1)*s_hy))/(2*Ly_top));
   xiy(1:ny,1)=2./((1 + 1j*gamma(2:ny+1)).*((1 + 1j*gamma(2:ny+1)) + (1 + 1j*gamma(1:ny))) );
   xiy(1:ny,2)=2./((1 + 1j*gamma(2:ny+1)).*((1 + 1j*gamma(2:ny+1)) + (1 + 1j*gamma(3:ny+2))) );
   
   % For Z Axis
   %-------------------------------------------------------------------
   gamma(1:npz_bot)           = cos((pi*(0:npz_bot-1)   *s_hz) /(2*Lz_bot));
   gamma(npz_bot+1:nz+2-npz_top) = 0;
   gamma(nz+3-npz_top:nz+2)   = cos((pi*(1-(nz+2-npz_top:nz+1)*s_hz))/(2*Lz_top));
   
   xiz(1:nz,1)=2./((1 + 1j*gamma(2:nz+1)).*((1 + 1j*gamma(2:nz+1)) + (1 + 1j*gamma(1:nz))) );
   xiz(1:nz,2)=2./((1 + 1j*gamma(2:nz+1)).*((1 + 1j*gamma(2:nz+1)) + (1 + 1j*gamma(3:nz+2))) );
   
end % create_pml
