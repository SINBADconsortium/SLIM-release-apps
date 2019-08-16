function [x,hst] = CRMN(A,b,x,opts)
%------------------------------------------------------------------------------
% Serial CRMN algorithm, supports block right hand sides
%
% USE:
%   [x,hst] = CRMN(A,b,x,opts)
%
% INPUT:
%   A                - hermitian matrix or spot operator
%   q                - right hand side, can be a matrix (each column is a separate rhs)
%   x                - initial guess
%   opts.tol         - CG tolerance, default = 1e-6
%   opts.maxit       - max. iterations, default = length(q)
%   opts.out_freq    - display output every # of iterations (default: 0, no output)
%
% OUTPUT:
%   x   - solution  
%   hst - history of residuals 
%   
% AUTHOR: Tristan van Leewen, 2014
%         Rafael Lago, 2014
%         Curt Da Silva, 2015
%         Mathias Louboutin, 2015
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
tol         = check_field(opts,'tol',1e-6);
maxit       = check_field(opts,'maxit',size(b,1));
out_freq    = check_field(opts,'out_freq',0);
nrhs        = size(b,2);

% Initialize CRMN pre-loop structures
%---------------------------------------------------
k  = 0;
if isempty(x)
    r = b; 
    x = zeros(size(b));
else
    r = b - A*x;
end

Ar = A*r;
p  = r;
Ap = Ar;

% 2-norm along the rows only
if nrhs==1
    innerprod = @(x,y) vec(x)'*vec(y);   
    n = @(x) norm(x);
else
    n = @(y) sum(abs(y).^2,1).^(1/2);
    innerprod = @(x,y) sum(conj(x) .* y,1);
end
nb  = n(b);
nr  = n(r);
nAp = n(Ap);
rAr    = innerprod(r,Ar);
hst    = max(nr./nb);

while true
    if nrhs==1
        gamma = rAr/nAp^2;
    else
        gamma = spdiags(vec(rAr./(nAp.^2)),0,nrhs,nrhs);    
    end
    
    x = x + p*gamma;
    r = r - Ap*gamma;
    nr = n(r);
    hstloc = max(nr./nb);
    
    hst = [hst hstloc];
    if out_freq > 0 && mod(k,out_freq)==0
        disp(['itr ' num2str(k,'%4.2d') ' | res ' num2str(hstloc,'%3.3e')]);
    end
    k = k+1;

    if (k>=maxit)||(hst(end)<=tol)
      break
    end
        
    Ar = A*r;
    oldrAr = rAr;
    
    rAr = innerprod(r,Ar);
    if nrhs==1
        beta = rAr/oldrAr;
    else
        beta = spdiags(vec(rAr./oldrAr),0,nrhs,nrhs);
    end
    p = r + p*beta;
    Ap = Ar + Ap*beta;
    nAp = n(Ap);   
end


end

