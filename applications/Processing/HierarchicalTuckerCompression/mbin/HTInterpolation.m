function [xsol, b, e,out] = HTInterpolation(fullData,dimTree,varargin)
%  HTINTERPOLATION - Convenience function for generating instances of tensor
%  interpolation problems + solving them with Hierarchical Tucker interpolation.
%  Subsamples the input data volume along specified dimensions,
%  then interpolates using HT interpolation and returns the result.
%
%  Curt Da Silva
%  HTOpt v0.1
%  curtd@math.ubc.ca
%
%  Usage:
%    [xsol, b, e, out] = HTInterpolation(fullData, dimTree, 'optionalArg1',optionalArgVal1, ...);
%
%  Input:
%    fullData           - full tensor to interpolate
%    dimTree            - corresponding dimension tree
%  
%  Optional Inputs:
%    'maxIter'          - maximum number of iterations (default: 100)
%    'tol'              - optimality tolerance (default: 1e-5)
%    'subsampling'      - percentage of points used to generate the training set (default: 0.1)
%    'subsample_dims'   - dimensions along which to subsample (default: 1:length(size(fullData)))
%    'rand_seed'        - random seed to use, for repeatability (default: 'default')
%    'skipOptimization' - if true, does not run the HT interpolation (useful for
%                       re-generating problem instances after the fact)
%    'logFile'          - file name to write the log file (default: [])
%    'method'           - one of 'SD', 'CG_PRP', 'CG', 'GN', see fitHT for more details
%    'noiseLevel'       - if > 0, adds random Gaussian noise at the specified level to
%                       the training data (default: 0)
%    'lambda'           - regularization value for the HT interpolation (default: 0)
%    'progTol'          - progress tolerance for the HT interpolation (default: 1e-6)
%    'verbosity'        - fitHT verbosity (default: 0 - silent)
%
%  Output:
%    xsol               - vectorized HT parameters of the estimated interpolant
%    b                  - training data
%    e                  - logical vector satisfying b(i) is known if and only if e(i)==1 for all i
%    out                - structure containing
%                         .solveTime - time to solve the interpolation problem
%                         .fhist     - vector of length #iterations,
%                                      fhist(k) = objective value at iteration k
%                         .trainErr  - relative training error
%                         .testErr   - relative test error
%                         .relErr    - relative (overall) error
    
    
   [maxIter,tol,subsampling,...
    subsample_dims,rand_seed,...
    skipOptimization,logFile,...
    method,noiseLevel,lambda,progTol,verbosity] = process_options(varargin, ...
                       'maxIter',100',...
                       'tol', 1e-4,...
                       'subsampling',0.1,...
                       'subsample_dims',1:length(size(fullData)), ...
                       'rand_seed','default', ...
                       'skipOptimization',false, ...
                       'logFile',[], ...
                       'method','GN', ...
                       'noiseLevel',0, ...
                       'lambda',0,...
                       'progTol',1e-6, ...
                       'verbosity',0);
      
   dims = size(fullData);
   if(~isempty(rand_seed))
       rng(rand_seed);
   end
   
   % Initial random point
   x0 = project(dimTree.randn(),dimTree);
   [U,B] = dimTree.fromVec(x0);
   B{1}{1} = B{1}{1}/norm(vec(B{1}{1}));
   x0 = dimTree.toVec(U,B);
  
   % Subsampling, removing points in the specified dimensions
   removeDims = prod(dims(subsample_dims));
   keepDims = prod(dims(setdiff(1:length(dims),subsample_dims)));
   numKeepShots = round(removeDims * subsampling);       
   J = randperm(removeDims, numKeepShots);
   
   R = opKron(opDirac(keepDims),opRestriction(removeDims,J));
   permutation = [subsample_dims,setdiff(1:length(dims),subsample_dims)];
   P = opPermute(dims,permutation);
       
   A = R* P;
   % Vector of 1s where there are known entries
   e = logical(A' * ones(size(A,1),1));

   b = fullData;
   b(~e) = 0;
   
   
   % Add noise if required
   if noiseLevel > 0
       n = sqrt(noiseLevel) * randn(prod(dims),1);
       n(~e) = 0;
       b = b + n;       
   end
   if skipOptimization
       xsol = []; out = []; return;
   end
   norm_b = norm(vec(b));
   
   %Interpolation
   out = struct;
   tic;
   [xsol,fk,fhist] = fitHT(logical(e),b,dimTree,...
                           'x0',x0,...
                           'maxIter',maxIter,...
                           'tol',tol,...
                           'method',method,...
                           'progTol',progTol,...
                           'verbosity',verbosity);
   out.solveTime = toc;
   out.fhist = fhist; 
   
   % Error output
   Y = dimTree.fullND(xsol); %Solution tensor
   Y(~e) = 0; 
   out.trainErr = norm(vec(Y) - vec(b))/norm(vec(b));
   
   Y = dimTree.fullND(xsol);
   Y(e) = 0;
   out.testErr = norm(vec(Y) - (vec(fullData - b)))/norm(vec(fullData - b));
   
   out.relErr = norm(dimTree.full(xsol) - vec(fullData));
                             

end
