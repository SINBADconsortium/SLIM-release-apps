%------------------------------------------------------------------------------
% Serial FGMRES algorithm 
%
% USE:
%   [x,hst] = FGMRES_MT(A,b,x,par)
%
% INPUT:
%   {A,idx}   - matrix in band storage format (no normalization needed)
%   b         - right hand side
%   x         - initial guess
%   par.tol   - FGMRES tolerance, default = 1e-6
%   par.maxit - For restarted GMRES: max. iterations, default = 10.
%   par.maxcy - Maximum number of cycles to perform. Default = inf.
%   par.precon- function handler to the preconditioner. It has to request only
%               one parameter. The call goes like w=precon(v). Default is
%               set to identity, that is, precon=@(x)(x).
%   flog      - file identifier for the (already open) log file. 
%               Let it unset or set it to 0 if you are targeting 
%               PERFORMANCE or simply do not want any message printed.
%
% OUTPUT:
%   x   - approximate solution  
%   hst - history of quasi-residuals 
% 
% AUTHOR: Rafael Lago (Principal)
%		  Mathias louboutin : Multiple RHS extension
%         Curt Da Silva : cleaned up code
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
function [x,hst] = FGMRES_MT(A,b,x,par)

% Default parameters
%--------------------------------------------------
ns=size(b,2);
tol = check_field(par,'tol',1e-6);
maxinnerit = check_field(par,'maxinnerit',10);
maxit = check_field(par,'maxit',inf);
n_threads = check_field(par,'n_threads',ns);
output_freq = check_field(par,'output_freq',0);

P = check_field(par,'precond',@(x) x);
if isa(P,'opSpot'), precond = @(x) P*x; else precond = P; end

%set_maxNumCompThreads(n_threads);
ompth=['OMP_SET_NUM_THREADS(' num2str(n_threads) ')'];
mtimesx('SPEEDOMP',ompth) ;
% Initialize stopping criterion
%--------------------------------------------------
% Compute initial residual
if ~isempty(x) && nnz(x)~=0
   r = b - A*x;
else
   r = b;
end

N = size(A,2);

k=1;
it_counter = 1;

blk_norm = @(x) sqrt(sum(abs(x).^2,1));
hst = 1; 
normr0 = blk_norm(r);
rnorm = normr0;
rnorm_hist = 1;
% Cycle loop
%--------------------------------------------------
while max(rnorm./normr0) > tol && floor(it_counter/maxinnerit)<=maxit
    %Initialization of Variables
    Z=zeros(N,maxinnerit,ns);
    V=zeros(N,maxinnerit+1,ns);
    H=zeros(maxinnerit+1,maxinnerit,ns);
    ej=zeros(maxinnerit+1,ns);
    w=zeros(N,1,ns);
    y=zeros(maxinnerit,1,ns);
    beta = blk_norm(r);
    ej(1,:) = 1;
    for i=1:ns, V(:,1,i) = r(:,i)/beta(i); end
   
    % Iteration loop
    %---------------
    while (k<=maxinnerit) && (hst(end)>tol)
        % Apply preconditioner
        Z(:,k,:)=precond(squeeze(V(:,k,:)));
        
        % Arnoldi
        %----------
        w(:,1,:)=A*squeeze(Z(:,k,:));
        % Classical Gram-Schmidt
        H(1:k,k,:)=mtimesx(V(:,1:k,:),'C',w);  
        w=w-mtimesx(V(:,1:k,:),H(1:k,k,:));
        
        % Modified Gram-Schmidt
        % for i=1:k
        %    H(i,k)=V(:,i)'*w;  
        %    w=w-H(i,k)*V(:,i);
        % end
        hstloc=0;
        for i=1:ns
            H(k+1,k,i)=norm(w(:,i));
            V(:,k+1,i)=w(:,i)/H(k+1,k,i);
            
            % Solve least squares and update solution  
            y(1:k,1,i)=H(1:k+1,1:k,i)\(beta(i)*ej(1:k+1,i)); 
            hstloc   =max(hstloc,norm(beta(i)*ej(1:k+1,i)-H(1:k+1,1:k,i)*y(1:k,1,i))/normr0(i));
        end
        
        hst = [hst hstloc];
        k = k + 1;
        it_counter = it_counter+1;
        
    end
    x = x+squeeze(mtimesx(Z(:,1:k-1,:),y(1:k-1,1,:))); 
	
    r = b - A*x;
    k=1;
    rnorm = blk_norm(r);
    rnorm_hist = [rnorm_hist, max(rnorm./normr0)];   
    if output_freq > 0 && mod(floor(it_counter/maxinnerit),output_freq)==0
        disp(['k ' num2str(floor(it_counter/maxinnerit),'%3.2d') ' res ' num2str(rnorm_hist(end),'%3.3e')]);
    end

end %Cycle loop

end
