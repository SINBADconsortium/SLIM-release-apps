function [x,hst] = CARPCG(R,idx,q,x,par)
%------------------------------------------------------------------------------
% Parallel CARPCG algorithm for band-storage matrices
%
% USE:
%   [x,hst] = CARPCG(R,idx,q,x,par)
%
% INPUT:
%   {R,idx} - normalized matrix in band storage format
%   q       - right hand side
%   x       - initial guess
%   par.w     - relaxation parameter, default = 1.5
%   par.tol   - CG tolerance, default = 1e-6
%   par.maxit - max. iterations, default = length(q)
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

spmd
   np = numlabs;
end
np = np{1};
if np==1
    error('Only one worker available, use serial version');
end

% Parse parameters
%--------------------------------------------------
w     = check_field(par,'w',1.5);
tol   = check_field(par,'tol',1e-6);
maxit = check_field(par,'maxit',length(q));
ns    = check_field(par,'ns',1);
n_threads = check_field(par,'n_threads',getenv('OMP_NUM_THREADS'));
if isempty(n_threads) ; n_threads = 1; end

n = par.size;
nxy = prod(n(1:end-1));

spmd
   %FIXME : REMOVE THIS!! It should always come in the same
   %format - the correct one.
   %--------------------------------------------------
   % Make sure diagonals of matrix are stored as rows
   %--------------------------------------------------
   if size(R,1)==length(idx)
      % do nothing
   elseif size(R,2)==length(idx)
      R = transpose(R);
   else
      error('Dimensions of R and idx do not match.');
   end
   
   % Distribute initial guess and preconditioned RHS
   %---------------------------------------------------
   x = pSPOT.pWindow.funWindowLast1HaloMakeCodist(x,n(end),np,1);
   q = pSPOT.pWindow.funWindowLast1HaloMakeCodist(q,n(end),np,1);

   b = pdsweep(R,idx,0*x,q,w,n,n_threads);
      assert(not(any(isnan(getLocalPart(b)))), ['CARPCG: an element of b, the'...
                           ' right hand side, is NaN after the initial sweep.']);

   % Initialize CARPCG pre-loop structures
   %---------------------------------------------------
   k = 0;
   r = b - x + pdsweep(R,idx,x,0*b,w,n,n_threads);
   p = r;

   normb = pnorm(b,nxy);
   normr = pnorm(r,nxy);
   
   % This variable will grow inside the loop, which is a priori inneficient...
   % but since its not used for computing and only for storing, I have the
   % feeling that it really doesn't matter. Please, someone correct me if I'm
   % wrong - Lago
   qres  = [normr];
   
   % Main loop
   %--------------------------------------------------
   while (k<maxit)&&(normr/normb>tol)    

      q = p - pdsweep(R,idx,p,0*b,w,n,n_threads);
      
      % alpha needs to be codistributed to make x codistributed also (as
      % opposed to a composite)
      alpha = codistributed(normr^2/pdot(p,q,nxy));
      
      x = x + alpha*p;
      r = r - alpha*q;
      
      normr = pnorm(r,nxy);
      beta  = codistributed(normr^2/qres(end)^2);
      
      p = r + beta*p;
      
      k = k + 1;
      qres = [qres normr];
      
   end
   
   % end of main loop
   %--------------------------------------------------
   x = pSPOT.pWindow.funWindowLast1HaloDropCodist(x,n(end),np,1);
end % spmd

hst = qres{1}/normb{1};

end % function CARPCG
