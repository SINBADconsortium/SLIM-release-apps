function [ x, fk, fhist ] = minFunc_hTuck( funObj, dimTree, varargin )
% MINFUNC_HTUCK - Minimizes a smooth function over the set of Hierarchical Tucker
% tensors, with optional regularization i.e.
%
%   min_{x} f(\phi(x)) + lambda * \sum_{t \in T} \|X^{(t)}\|_{F}^2
%                                              + \|(X^{(t)})^{\dagger}\|_{F}^2
%
% See 'Optimization on the Hierarchical Tucker manifold -
% applications to tensor completion', C. Da Silva and F. Herrmann,
% 2013, for more details. 
%
% Curt Da Silva
% HTOpt v 1.0
% curtd@math.ubc.ca
%
%
% Usage:
%   [x, fk, fhist] = minFunc_hTuck( funObj, dimTree, 'optionalInput1',optionalInput1,... )
% 
% Input:
%   funObj      - function handle which computes the objective and
%                 the Euclidean gradient at specified point on the manifold
%   dimTree     - dimension tree object specfying the HT format
% 
% Optional Input:
%   maxIter     - maximum number of iterations (default: 100)
%   x0          - initial guess for x0 (default: random)
%   tol         - optimality tolerance (default: 1e-5)
%   verbosity   - 0 : silent mode, no output
%               - 1 : per-iteration output (default)
%               - 2 : verbose output
%   logFile     - file name to save output (default: [])
%   method      - 'SD' : steepest descent (very slow, NOT recommended)
%               - 'CG' : CG-Descent
%               - 'CG_PRP' : Polyak Ribiere Conjugate Gradient
%               - 'GN' : Gauss-Newton method (default)
%   progTol     - Tolerance to ensure optimization progress (default: 1e-9)
%   suffDec     - Armijo sufficient descent parameter (default: 0.1)
%   theta       - Step size decrease parameter (default: 0.5)
%   maxLS       - Maximum number of line search iterations (default: 50)
%   lambda      - Regularization parameter (default: 0)
%   distributed - if true, full HT tensors will be distributed
%                 vectors (default: false)   
%   dimtreeFunc - if true, objective accepts HT parameters x,
%                 parameters, outputs objective, Riemannian gradient
%                 (default: false)
%   outputUpdate - print output every # of iterations (default: 1)
%
%  Output:
%     x         - solution HT parameters
%     fk        - final objective value (including regularization,
%                 if any)
%     fhist     - array of per-iteration objective values 

VERBOSE = 2;
ITER = 1;
SILENT = 0;

[maxIter, x0, tol, verbosity, logFile, method, ...
    progTol, delta, theta, maxLS,lambda, ... 
    distributedMode,dimtreeFunc,outputUpdate] =  ...
    process_options(varargin, ...
    'maxIter',100, ...
    'x0',[],...
    'tol',1e-5, ...
    'verbosity', ITER, ...
    'logFile',[],...
    'method','GN',...
    'progTol',1e-5, ...
    'suffDec',1e-2, ... 
    'theta',0.5, ...
    'maxLS',50,...
    'lambda',0,...
    'distributed',false,...
    'dimtreeFunc',false,...
    'outputUpdate',1 );

METHODS = {'SD','CG','CG_PRP','GN'};
SD = 1;
CG = 2;
CG_PR = 3;
GN = 4;
METHOD = find(cellfun(@(s) strcmp(s,upper(method)),METHODS));
if isempty(METHOD)
    error('Must provide a valid method');
end

if ~isempty(logFile)
    log_fid = fopen(logFile,'w');
else
    log_fid = -1;
end

    function output(str, minVerbosity)
        if(verbosity >= minVerbosity)
            disp(str);
            if(log_fid > 0)
                fprintf(log_fid,strcat(str,'\n'));               
            end
        end
    end

if isempty(x0)        
    x = project(dimTree.randn(),dimTree);
    [U,B] = dimTree.fromVec(x);
    B{1}{1} = B{1}{1}/norm(vec(B{1}{1}));
    x = dimTree.toVec(U,B);
else
    x = project(x0,dimTree);
end

    function [fk,gk] = htuckObjective(x,computeGradient, fk,gk)
    % Encapsulates objective/Riemannian gradient evaluation + regularization on the
    % HT manifold 
        if exist('computeGradient','var')==0
           computeGradient = false; 
        end
        
        if dimtreeFunc
            if computeGradient
                [fk,gk] = funObj(x);
                gevals = gevals + 1;
            else
                fk = funObj(x);                
            end
            fevals = fevals+1;
        else        
            if exist('fk','var')==0                       
                if distributedMode
                    [fk,gk] = funObj(dimTree.fullDist(x));
                else
                    [fk,gk] = funObj(dimTree.full(x));
                end
                fevals = fevals + 1;                  
            end
            
            if ~isempty(lambda) & lambda > 2*eps
                if computeGradient
                    [f2,g2] = gramian_regularizer(dimTree,x,lambda);  
                else
                    f2 = gramian_regularizer(dimTree,x,lambda);
                end
                fk = fk + f2;
            end        
            
            if computeGradient
                if distributedMode
                    J = oppHTuckJ2(dimTree,x);
                else
                    J = opHTuckJ2(dimTree,x);
                end
                gk = J' * gk; gevals = gevals + 1; 
                if ~isempty(lambda) & lambda > 2*eps
                    gk = gk + g2;
                end
            end      
        end
    end
    
fevals = 0;
gevals = 0;

% Compute objective + full HT gradient
[fk,gk] = htuckObjective(x,true);

alpha = 1;
fhist = zeros(maxIter,1);
fhist(1) = fk;
resetStep = false;

for itr=1:maxIter
%%--------- Search direction    
    innerprod = @(dx,dy) vec(dx)' * vec(dy);
    gk_norm = norm(gk);
    if itr ==1 || METHOD == SD
        pk = -gk; Lk = 1;
    else
        if resetStep 
            pk = -gk; Lk = 1; resetStep = false;
        else
            yk = gk - gkprev;
            sk_norm = norm(sk);
            yk_norm = norm(yk);
            yk_sk = innerprod(yk,sk);
            pkprev = sk/alpha;
            Lk = real(yk_sk)/sk_norm^2;
            switch METHOD
              case GN
                HGN = opHTuckGN(dimTree,x);
                pk =  HGN \ (-gk);
              case CG
                %  CG scheme 1.2/1.3/1.5/1.6 in 'A New Conjugate Gradient
                %Method with Guaranteed Descent and an Efficient Line
                %Search'
                etak = -1 /(sk_norm * min(0.01,norm(gkprev)));
                beta = real(innerprod(yk - 2*sk * yk_norm^2 / yk_sk, gk)) / real(yk_sk);
                beta = max(beta,etak);
                pk = -gk + beta * sk;
              case CG_PR
                %Hestenes-Stiefel CG           
                beta = real(innerprod(yk,gk))/real(yk_sk);   
                pk = -gk + beta * sk;             
            end
        end
    end
    gk_pk = real(innerprod(gk,pk));
    if gk_pk > -1e-9
       output('Uphill direction, resetting',VERBOSE);
       pk = -gk; gk_pk = -gk_norm^2;
    end
    
%%--------- Output
    out_string = ['k ' num2str(itr,3) ' fk: ' num2str(fk,'%3.5e') ...
        ' g : ' num2str(gk_norm,'%3.3e') ...
        ' optim. cond: ' num2str(abs(gk_pk),'%3.3e') ...
        ' fevals: ' num2str(fevals) ' alpha ' num2str(alpha,'%3.3e') ];
    if mod(itr-1,outputUpdate)==0
        output(out_string, ITER);
    end
%%--------- Termination check
    if itr > 2 &&  abs(fk - fkprev) < progTol * abs(fk)
        if alpha < 1e-6
            pk = -gk; gk_pk = -gk_norm^2;
            output('Direction reset',ITER);
        else
            output('Function value changing by less than tolerance',ITER);
            break;
        end
    end
    
    if itr > 2 && (gk_norm < tol || sqrt(abs(alpha)) < tol || (alpha > 1e-3 && alpha * abs(gk_pk) <= 1e-20 * fk))
        output('Directional derivative/gradient less than tolerance',ITER);
        break;
    end
    
%%--------- Line search    
    if Lk > 0 
        alpha = abs(gk_pk)/(Lk * norm(pk)^2);        
    end
    if alpha < 1e-6 || alpha > 1e6 
        alpha = 1;
    end
    fmin = inf;
    alphamin = alpha;
    acceptedStep = false;
    for i=1:maxLS
        xnext = x + alpha * pk;
        if i > 1
            fnextprev = fnext;
        end
        fnext = htuckObjective(xnext);
        output(['LS alpha ' num2str(alpha) ' phi(alpha) ' num2str(fnext,'%3.6e')],VERBOSE);
        
        fmin = min(fmin,fnext);
        if fmin == fnext
            alphamin = alpha;
        end        
        if i==1 && (fnext - fk < delta * gk_pk * alpha)
            %Initial step works: look for a larger step size that works
            output(['Initial step is Armijo, increasing step size'],VERBOSE);
            fnext_init = fnext;
            
            step_size_init = alpha;
            for j=2:maxLS
                alpha = alpha / theta;
                xnext = x + alpha*pk;
                fnextprev = fnext;
                fnext = htuckObjective(xnext);
                output(['LS alpha ' num2str(alpha) ' phi(alpha) ' num2str(fnext,'%3.6e')],VERBOSE);
                
                if fnext - fk >= delta * gk_pk * alpha || fnextprev < fnext
                    break;
                end
                fmin = min(fmin,fnext);
                if fmin == fnext
                    alphamin = alpha;
                end
                
            end
            if j >= 3 %Backtracking past the initial point is just a waste now
                break;
            end
            alpha = step_size_init;
            fnext = fnext_init;
            
        end
        if i > 1 && (fnextprev < fnext || abs(fnextprev - fnext) < abs(fnext) * 1e-8)
            %            resetStep = true;           
            %break; %stagnation
        end
        fmin == min(fmin,fnext);
        if fmin == fnext
            alphamin = alpha;
        end
        
        if acceptedStep
            if fnext - fk >= delta * gk_pk * alpha || ...
                    fnextprev < fnext                
                break; %Forwardtracked too much
            end
        else
            if fnext - fk < delta * gk_pk * alpha
                acceptedStep = true; %Armijo, see if the next step is still Armijo
            end
        end
        alpha = alpha * theta;
    end
    
%%--------- Update variables + vector transport
    alpha = alphamin;
    x = project(x + alpha * pk,dimTree);
    fkprev = fk;
    gkprev = project_horizontal(x,gk,dimTree);
    sk = project_horizontal(x,alpha*pk,dimTree);       
    
%%--------- Function + gradient update
    [fk,gk] = htuckObjective(x,true);    
    
    fhist(itr+1) = fk;
end
fhist = fhist(1:itr);

end

