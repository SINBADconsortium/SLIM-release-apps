%------------------------------------------------------------------------------
% Serial FGMRES algorithm for band-storage matrices
%
% USE:
%   [x,res,hst] = FGMRES(A,b,x,par)
%
% INPUT:
%   A         - matrix or SPOT operator that can perform matrix-vector products
%   b         - right hand side
%   x         - initial guess
%   par.tol   - FGMRES tolerance, default = 1e-6
%   par.maxinnerit - For restarted GMRES: max. iterations, default = 10.
%   par.maxit - Maximum number of cycles to perform. Default = inf.
%   par.precon- function handler to the preconditioner. It has to request only
%               one parameter. The call goes like w=precon(v). Default is
%               set to identity, that is, precon=@(x)(x).
%   par.output_freq - display output every # of outer iterations (default: 0, no output)
%
% OUTPUT:
%   x   - approximate solution  
%   res - history of residuals
%   hst - history of quasi-residuals 
% 
% AUTHOR: Rafael Lago
%         Igor Fonseca
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%      
% Date: August, 2014
%
% Updated : Curt Da Silva, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%------------------------------------------------------------------------------
function [x,res,hst] = FGMRES(A,b,x,par)

if size(x,2) > 1
    [x,res,hst] = FGMRES_MT(A,b,x,par); return; 
end

% Parse parameters
%--------------------------------------------------

tol = check_field(par,'tol',1e-6);
maxinnerit = check_field(par,'maxinnerit',10);
maxit = check_field(par,'maxit',inf);
output_freq = check_field(par,'output_freq',0);

P = check_field(par,'precond',@(x) x);
precond_dirac = false;
if isa(P,'opSpot')
    if isa(P,'opDirac'), precond_dirac = true;
    else precond = @(x) P*x; end
else precond = P; 
end
% Initialize stopping criterion
%--------------------------------------------------
% Compute initial residual
if nnz(x)~=0
   r = b - A*x;
else
   r = b;
   x = zeros(size(b));
end

N = size(A,2);

it_counter = 1;
normr0     = norm(r);
hst        = 1; % History of quasi-residuals - first is always one
ej=zeros(maxinnerit+1,1);
ej(1) = 1;   
res = [];
% Cycle loop
%--------------------------------------------------
while hst(end)>tol && floor(it_counter/maxinnerit)<maxit
    k=1;
    %Initialization of Variables
    if precond_dirac
        %clear Z;
        %Z(N,maxinnerit+1) = 0;
        Z = zeros(N,maxinnerit);
    else
        %clear Z V;
        %Z(N,maxinnerit) = 0;
        %V(N,maxinnerit+1) = 0;
        Z = zeros(N,maxinnerit);
        V = zeros(N,maxinnerit+1);
    end
    
    H=zeros(maxinnerit+1,maxinnerit);

    y=zeros(maxinnerit,1);    
    beta = norm(r);
    if precond_dirac
        Z(:,1)=r/beta;
    else
        V(:,1)=r/beta;
    end    
        
    % Iteration loop
    %---------------
    while (k<=maxinnerit) && (hst(end)>tol)
        
        % Apply preconditioner
        if ~precond_dirac
            Z(:,k)=precond(V(:,k));
        end
        % Arnoldi
        %----------
        w=A*Z(:,k);

        % Classical Gram-Schmidt
        if precond_dirac
            H(1:k,k)=Z(:,1:k)'*w;
            w = w-Z(:,1:k)*H(1:k,k);
        else
            H(1:k,k)=V(:,1:k)'*w;  
            w=w-V(:,1:k)*H(1:k,k);
        end
        
        H(k+1,k)=norm(w);
        if precond_dirac
            Z(:,k+1) = w/H(k+1,k);
        else
            V(:,k+1) = w/H(k+1,k);
        end
        % Solve least squares and update solution  
        y(1:k,1)=H(1:k+1,1:k)\(beta*ej(1:k+1,1)); 
        
        hst   = [hst norm(beta*ej(1:k+1,1)-H(1:k+1,1:k)*y(1:k,1))/normr0];
        k = k + 1;
        it_counter = it_counter+1;        
    end
   
    x = x+Z(:,1:k-1)*y(1:k-1,1); 
    r = b - A*x;
    if output_freq > 0 && mod(ceil(it_counter/maxinnerit),output_freq)==0
        disp(['k ' num2str(floor(it_counter/maxinnerit),'%3.2d') ' res ' num2str(norm(r)/normr0,'%3.3e')]);
    end    
    if norm(r)/normr0 < tol, break; end
    res = [res, norm(r)]; 
end %Cycle loop
end  

