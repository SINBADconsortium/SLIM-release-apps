function [ U, B ] = fitTucker( trainIndices, rhs, dims, ranks, varargin )
%FITTUCKER - Fits a tensor with missing entries in the Tucker
%format. Implements the method from 'Low rank tensor completion by 
%Riemannian optimization', D. Kressner, M. Steinlechner,
%B. Vandereycken, 2013.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
%
% Usage:
%   [U,B] = fitTucker(trainIndices, rhs, dims, ranks, 'optionalArg1','optionalArgVal1',...);
%   X_recovered = ttm(B,U,1:length(U));
% 
% Input:  
%    trainIndices - 1/0 logical vector of length == length(rhs), which should satisfy
%                   trainIndices(i) == true if and only if rhs(i) is known
%    rhs          - zero-filled training data
%    dims         - dimensions of the target tensor
%    ranks        - Tucker ranks of the target tensor (need length(ranks) == length(dims) )
%    
% Optional inputs:
%    'maxIter'    - maximum number of iterations (default: 100)
%    'x0'         - initial guess for optimization (default: x0 is random)
%    'tol'        - optimality tolerance (default: 1e-4)
%    'verbosity'  - 0: no iteration output
%                 - 1: per-iteration output (default)
%                 - 2: verbose output
%    'method'     - 'SD' - Steepest Descent 
%                 - 'CG_PR' - Polyak Ribiere Conjugate Gradient (default)
%    'theta'      - step size decrease parameter (default: 0.5)
%    'maxLS'      - maximum number of line search iterations (default: 50)
%    'randseed'   - random seed - to ensure repeatability (default: 'default')
%
% Output:
%    [U,B]        - recovered Tucker parameters, which can be expanded to a full tensor via ttm(B,U,1:length(U))
SILENT = 0;
ITER = 1;
VERBOSE = 2;
    
[maxIter, x0, tol, verbosity,  method, ...
    progTol, delta, theta, maxLS] =  ...
    process_options(varargin, ...
    'maxIter',100, ...
    'x0',[],...
    'tol',1e-4, ...
    'verbosity', 1, ...
    'method','CG_PR',...
    'progTol',1e-6, ...
    'suffDec',1e-4, ... 
    'theta',0.5, ...
    'maxLS',20);

SD = 1;
CG_PR = 2;
METHODS = {'SD','CG_PR'};
METHOD = find(cellfun(@(s) strcmp(s,upper(method)),METHODS));

zeroIndices = ~trainIndices;
rhs_norm = norm(vec(rhs));
function v = vecParams(U,B)
% Vectorizes the cell array parameters in to a vector
    v = [];
    for j=1:length(U)
        v = [v;vec(U{j})];
    end
    v = [v;vec(B)];
end
function [U,B] = unvecParams(x)
% Reshapes vectorized paramters in to cell arrays
    numel = 0;
    U = cell(d,1);
    for j=1:d
        ni = dims(j); ri = ranks(j);
        U{j} = reshape(x(numel+1:numel+ni*ri),[ni,ri]);
        numel = numel + ni*ri;
    end
    B = reshape(x(numel+1:end),ranks);
    
end

function output(str,verb)
    if verb <= verbosity
        disp(str);
    end
end
    
    
d = length(dims);
rhs = reshape(rhs,dims);

% Initialization
if isempty(x0)
    U = cell(d,1);
    B = randn(ranks);
    for i=1:d
        [U{i},R] = qr(randn( dims(i),ranks(i) ),0);
        B = ttm(B,R,i);
    end   
    B = B ./ norm(vec(B));
    x = vecParams(U,B);
else
    x = x0;
    [U,B] = unvecParams(x);
end
      
function X = fullTensor(x)
% Expands the tensor parameters in to a full tensor
    [U,B] = unvecParams(x);
    X = ttm(B,U,1:d);
end
function dX = fullGradient(x,dx)
% Expands the Riemannian gradient in to a full tensor (for the line search)
    [U,B] = unvecParams(x);
    [dU,dB] = unvecParams(dx);
    dX = ttm(dB,U,1:d);
    for j=1:d
        Ut = U;
        Ut{j} = dU{j};
           Z = ttm(B,Ut,1:d);
           dX = dX + Z; 
    end       
end
function [f,g,r] = obj(x)  
% Objective function + Riemannian gradient
    X = fullTensor(x);
    [U,B] = unvecParams(x);
    X(zeroIndices) = 0;
    r = X - rhs;
    f = 0.5 * norm(vec(r))^2;
    if nargout > 1
        [dU,dB] = project_tucktangent(U,B,r); 
        g = vecParams(dU,dB);
    end
end
function y = xplusay(x,dx,alpha)
% HOSVD-based retraction exploiting the Tucker structure of the problem
    [U,B] = unvecParams(x);
    [dU,dB] = unvecParams(dx);
    S = zeros(2*ranks);
    r = ranks;
    idx = [''];
    for j=1:d
        idx = [idx,'1:' num2str(r(j)) ','];
    end
    idx = idx(1:end-1);
    eval(['S(' idx ') = B + alpha * dB;']);
    
    for j=1:d
        idx = [''];
        for k=1:d
            if j~=k
                idx = [idx,'1:' num2str(r(k)) ','];
            else
                idx = [idx, num2str(r(k)+1) ':end,'];
            end
            end
            idx = idx(1:end-1);
            eval(['S(' idx ') = alpha * B;']);
    end
    Q = cell(length(U),1);
    for j=1:d
        [Qi,R] = qr([U{j} dU{j}]); 
        Q{j} = Qi;
        
        S = ttm(S,R,j);
    end
    [V,C] = HOSVD(S,r);
    for j=1:d           
        Q{j} = Q{j} * V{j};                       
    end
    y = vecParams(Q,C);
end

x = vecParams(U,B);
xi = 0;
[fk,gk,r] = obj(x);

stepSize = 1;

function z = innerprod(x,dx,dy)
% Inner product between two horizontal vectors
    [U,B] = unvecParams(x);
    [dU,dB] = unvecParams(dx);
    [dV,dC] = unvecParams(dy);
    z = vec(dB)' * vec(dC);
    for j=1:d
        z = z + vec(B)' * vec(ttm(B,dV{j}' * dU{j},j));
    end
end

X = [];
t = [];
Xold = [];
value_old = [];
for itr=1:maxIter
    norm_gk = sqrt(innerprod(x,gk,gk));
    out_string = ['k ' num2str(itr) ' f ' num2str(fk,'%3.6e') ' g ' num2str(norm_gk,'%3.3e') ' alpha ' num2str(stepSize,'%3.3e')];
    output(out_string,ITER);
            
    % Search direction 
    if itr == 1 || METHOD == SD
        pk = -gk;  
        gk_pk = -norm_gk^2;
    else        
        
        theta = innerprod(x,gkprev,gk)/innerprod(x,gk,gk);
        if(theta >= 0.1)
            pk = -gk;
            gk_pk = -norm_gk^2;
        else        
            beta = max(innerprod(x,gk,gk - gkprev)/innerprod(x,gkprev,gkprev),0);  
            output(['CG beta ' num2str(beta,'%3.3e')],VERBOSE);
            pk = -gk + beta * pkprev;

            gk_pk = innerprod(x,gk,pk);
              
        end
    end
    
    % Converged
    if norm_gk < tol         
        break;
    end
    if itr > 1 && abs(sqrt(fk) - sqrt(fkprev)) < progTol * sqrt(fk);
        output('Insufficient function decrease',ITER); 
        break;
    end
    
    
    %Linearized approximate minimizer for least squares
    dX = fullGradient(x,pk);
    dX(zeroIndices) = 0;
    stepSize = max(-(vec(dX)' * vec(r))/norm(vec(dX))^2, 1e-10);
    if itr > 1    
        for i=1:maxLS
            xnew = xplusay(x,pk,stepSize);        
            fknew = obj(xnew);
            output(['alpha ' num2str(stepSize,'%3.3e') ' phi(alpha) ' num2str(fknew,'%3.6e')],VERBOSE);
            if fknew - fk < delta * gk_pk * stepSize           
               break; 
            end
            stepSize = stepSize * theta;
        end
    else
        xnew = xplusay(x,pk,stepSize);
    end
    [dUold,dBold] = unvecParams(gk);
    [Uold,Bold] = unvecParams(x);
    [U,B] = unvecParams(xnew);
    [dU,dB] = transport_tucktangent(Uold,Bold,dUold,dBold,U,B);
    gkprev = vecParams(dU,dB);
    [dUold,dBold] = unvecParams(pk);
    [dU,dB] = transport_tucktangent(Uold,Bold,dUold,dBold,U,B);
    pkprev = vecParams(dU,dB);
    
    x = xnew;
    fkprev = fk;
    [fk,gk,r] = obj(x);
end

end


function [dU,dB] = project_tucktangent(U,B,X)
% Extracts the Riemannian gradient at (U,B) from a tangent vector X
    ranks = size(B);
    dU = cell(length(U),1);
    dB = X;
    d = length(U);
    for i=1:d
       dB = ttm(dB,opMatrix(U{i}' ),i);
    end
    for i=1:d
       Binv = pinv(matricize(B,i,setdiff(1:d,i)));       
       Z = opProjection(U{i},true,true);         
        
       Ut = U;
       for j=1:d
          if i==j
             Ut{j} = opDirac(size(Ut{j},1));
          else
             Ut{j} = Ut{j}';
          end
       end
       
       Q = ttm(X,Ut,1:d);
       Q = Z * matricize(Q,i) * Binv;
              
       dU{i} = Q;
       
    end
end

function [dU,dB] = transport_tucktangent(Uold,Bold,dUold,dBold,U,B)
% Transports the tangent vector (dUold,dBold) from the horizontal space at (Uold,Bold) to 
% the horizontal space at (U,B)
    d = length(U);
    UtdV = cell(length(dUold),1);
    sz = zeros(length(d),1);
    for i=1:d
        UtdV{i} = U{i}' * Uold{i};
        sz(i) = size(U{i},1);
    end
    dB = ttm(dBold,UtdV,1:d);
    
    for i=1:d
        Ut = U;
        for j=1:d
            if i==j
                Ut{j} = U{j}' * dUold{j};
            else
                Ut{j} = U{j}' * Uold{j}; 
            end
        end
        dB = dB + ttm(Bold,Ut,1:d);
    end
    
    dU = cell(length(dUold),1);
    for i=1:d
       Z = opProjection(U{i},true,true); 
       
       Ut = U;
       for j=1:d
          if i~=j                        
             Ut{j} = U{j}' * Uold{j};
          end
       end
       Y = ttm(dBold,Ut,1:d);
              
       for j=1:d
          Ut = U;
          for k=1:d
              if k==j
                  if i==j
                      Ut{k} = dUold{k};
                  else
                      Ut{k} = U{k}' * dUold{k};
                  end                 
              else
                  if i==k
                      Ut{k} = Uold{k};
                  else
                      Ut{k} = U{k}' * Uold{k};
                  end                  
              end
          end
          Y = Y + ttm(Bold,Ut,1:d);          
       end
       
       Binv = pinv(matricize(B,i,setdiff(1:d,i)));
       dU{i} = Z * matricize(Y,i,setdiff(1:d,i)) * Binv;              
    end            
    
end

