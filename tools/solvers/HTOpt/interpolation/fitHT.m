function [ x, fk, fhist ] = fitHT( trainIndices, rhs, dimTree, varargin )
% FITHT - Fits a set of (missing) tensorial data using the Hierarchical Tucker
% format corresponding to the provided dimension tree, i.e.
%
%  min_{x} 0.5 * \|P_{\Omega} \phi(x) - b\|_2^2 + lambda * \sum_{t \in T} \|X^{(t)}\|_{F}^2
%                                              + \|(X^{(t)})^{\dagger}\|_{F}^2
%
% See 'Optimization on the Hierarchical Tucker manifold -
% applications to tensor completion', C. Da Silva and F. Herrmann,
% 2013, for more details. 
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   [x, fk, fhist] = fitHT( trainIndices, rhs, dimTree, 'optionalInput1',optionalInput1,... )
% 
% Input:
%   trainIndices - 1-0 logical vector with length == length(rhs). Ones
%                 correspond to the locations where there is data
%                 in rhs.
%   rhs         - input data (vector) with zeros
%   dimTree     - dimension tree object 
% 
% Optional Input:
%   maxIter     - maximum number of iterations (default: 100)
%   x0          - initial guess for x0 (default: random)
%   tol         - optimality tolerance (default: 1e-5)
%   verbosity   - 0 : silent mode, no output
%               - 1 : per-iteration output (default)
%               - 2 : verbose output
%   logFile     - file name to save output (default: [])
%   method      - 'SD' : steepest descent (NOT recommended)
%               - 'CG' : Conjugate Gradient
%               - 'CG_PRP' : Polyak Ribiere Conjugate Gradient
%               - 'GN' : Gauss-Newton method (default)
%   progTol     - Tolerance to ensure optimization progress (default: 1e-9)
%   suffDec     - Armijo sufficient descent parameter (default: 0.1)
%   theta       - Step size decrease parameter (default: 0.5)
%   maxLS       - Maximum number of line search iterations (default: 50)
%   lambda      - Regularization parameter (default: 0)

%  Output:
%     x         - solution HT parameters
%     fk        - final objective value
%     fhist     - array of per-iteration function values
if size(rhs,2) > 1 || length(size(rhs)) > 2
    rhs = vec(rhs);
end
    
zeroIndices = ~trainIndices;

funObj = @(x) LSLinearMisfit(x,rhs,zeroIndices);

args = varargin;
args{end+1} = 'distributed'; args{end+1} = isdistributed(rhs) || iscodistributed(rhs);

[x,fk,fhist] = minFunc_hTuck( funObj, dimTree, args{:} );


end


