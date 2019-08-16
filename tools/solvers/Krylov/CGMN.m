function [x,hst] = CGMN(A,b,x,par)
%------------------------------------------------------------------------------
% Serial CGMN algorithm for band-storage matrices
%
% USE:
%   [x,hst] = CGMN(A,b,x,par)
%
% INPUT:
%   A - normalized matrix in band storage format
%   b       - right hand side
%   x       - initial guess
%   par.w     - relaxation parameter, default = 1.5
%   par.tol   - CG tolerance, default = 1e-6
%   par.maxit - max. iterations, default = length(b)
%   par.ns    - number of double sweeps per iterations, default = 1
%
% OUTPUT:
%   x   - solution  
%   hst - history of residuals 
%   
% AUTHOR: Tristan van Leewen
%         Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: May, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------

% Parse parameters
%--------------------------------------------------
tol   = check_field(par,'tol',1e-6);
maxit = check_field(par,'maxit',length(b));
out_freq = check_field(par,'out_freq',0);

nrhs  = size(b,2);
% Construct preconditioned RHS and sweeping operator
%---------------------------------------------------

% Initialize CGMN pre-loop structures
%---------------------------------------------------
k = 0;
r = b - A*x;
p = r;

norm = @(y) sum(abs(y).^2,1).^(1/2);

normb = norm(b);
normr = norm(r);


% This variable will grow inside the loop, which is a priori inneficient...
% but since its not used for computing and only for storing, I have the
% feeling that it really doesn't matter. Please, someone correct me if I'm
% wrong - Lago
hst   = [normr/normb];

% Main loop
%--------------------------------------------------
while true
   Ap    = A*p;
   pAp   = sum(conj(p).*Ap,1);
   alpha = spdiags(vec(normr.^2./pAp),0,nrhs,nrhs);
   
   x     = x + p*alpha;
   r     = r - Ap*alpha;

   normr_old = normr;
   normr     = norm(r);
   hst   = [hst max(normr./normb)];
   k = k + 1;
   
   if (k>=maxit)||(hst(end)<=tol)
      break
   end
   
   beta  = spdiags(vec(normr./normr_old).^2,0,nrhs,nrhs);   
   
   p     = r + p*beta;
   if out_freq > 0 && mod(k-1,out_freq)==0
       disp(['k ' num2str(k) ' res ' num2str( max(normr./normb),'%3.3e')]);
   end

end
% end of main loop
%--------------------------------------------------

end


