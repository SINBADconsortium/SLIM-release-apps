function [x,hst] = CRMN_old(R,idx,q,x,par)
%------------------------------------------------------------------------------
% Serial CRMN algorithm for band-storage matrices
%
% use:
%   [x,hst] = CRMN(R,idx,q,x,par)
%
% input:
%   {R,idx} - normalized matrix in band storage format
%   q       - right hand side
%   x       - initial guess
%   par.w     - relaxation parameter, default = 1.5
%   par.tol   - CR tolerance, default = 1e-6
%   par.maxit - max. iterations, default = length(q)
%   par.ns    - number of double sweeps per iterations, default = 1
%   flog      - file identifier for the (already open) log file. 
%               Let it unset or set it to 0 if you are targeting 
%               PERFORMANCE or simply do not want any message printed.
%
% output:
%   x   - solution  
%   hst - history of residuals  
%------------------------------------------------------------------------------  

% Parse parameters
%--------------------------------------------------
w     = check_field(par,'w',1.5);
tol   = check_field(par,'tol',1e-6);
maxit = check_field(par,'maxit',length(q));
ns    = check_field(par,'ns',1);
n_threads = check_field(par,'n_threads',getenv('OMP_NUM_THREADS'));
if isempty(n_threads) ; n_threads = 1; end


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
%--------------------------------------------------

% Construct preconditioned RHS and sweeping operator
%---------------------------------------------------
b = dsweep_old(R,idx,w,0*x,q,ns,n_threads);
Afun = @(x)(x - dsweep_old(R,idx,w,x,0*b,ns,n_threads));

% Initialize CRMN pre-loop structures
%---------------------------------------------------
k=0;
r = b - Afun(x);
Ar = Afun(r);
p = r;
Ap = Ar;
rAr = (r'*Ar);

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
   gamma = rAr/(norm(Ap)^2);

   x = x + gamma*p;
   r = r - gamma*Ap;

   normr = norm(r);
   hst = [hst normr/normb];
   k = k + 1;
   
   if (k>=maxit)||(hst(end)<=tol)
      break
   end
   
   Ar = Afun(r);
   
   oldrAr = rAr;
   rAr = (r'*Ar);
   beta = rAr/oldrAr;

   p = r + beta*p;
   Ap = Ar + beta*Ap;   
end
% end of main loop
%--------------------------------------------------

end