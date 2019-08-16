%------------------------------------------------------------------------------
% Serial Jacobi relaxation for band-storage matrices
%
% USE:
%   x = Jacobi_relax(A,idx,b,x,par)
%
% INPUT:
%   {A,idx}   - matrix in band storage format (no normalization needed)
%   b         - right hand side
%   x         - initial guess
%   par.tol   - FGMRES tolerance, default = 1e-6
%   par.maxit - Maximum number of iterations to be performed. Usually this is
%               set to one or two
%   par.omega - relaxation parameter
%
% OUTPUT:
%   x   - approximate solution  
% 
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: August, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function x = Jacobi_relax(A,idx,b,x,par)

% The standard Jacobi relaxation is usually writen as
%    xₖ₊₁ = ωD⁻¹(b - Rxₖ) + (1-ω)xₖ
% with A = R+D. However, noting that R = A - D, we find that
%    xₖ₊₁ = ωD⁻¹(b - Axₖ + Dxₖ) + (1-ω)xₖ = ωD⁻¹rₖ + ωxₖ + (1-ω)xₖ = ωD⁻¹rₖ + xₖ
% which is the formula we use down here. Notice that D=transpose(A(par.NNN,:))
% We do the first iteration outside the loop, because we want to avoid the 
% residual computation if possible.

nt = size(x,1);
Dinv = @(u)(u./reshape(A(par.NNN,:),nt,1));
if nnz(x)~=0
   r = b - Hmvp(A,idx,x);
else
   r = b;
end
x = par.omega*Dinv(r) + x;
% end of first iteration

for i=1:par.maxit-1
   r = b - Hmvp(A,idx,x);
   x = par.omega*Dinv(r) + x;
end

end
