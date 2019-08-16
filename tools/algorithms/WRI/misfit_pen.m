function [f,g,h] = misfit_pen(m,Q,D,model,params)

% penalty misfit for FWI
%
% use:
%   [f,g,h] = misfit_pen(m,Q,D,model,params)
%
% input:
%   m  - model parameters (slowness squared)
%   Q  - source matrix
%   D  - data
%   model - struct with model paramaters
%   params.lambda - tradeoff parameter (scalar) between PDE and data misfit
%    
%
% output:
%   f - function value
%   g - gradient
%   h - diagonal approximation of Hessian
%   f_aux - auxilary misfit measures: separate PDE and data misfit in a
%   structure
%

% Author: Tristan van Leeuwen, Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: April, 2014
%
% Updated by: Curt Da Silva, 2015/2016
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if exist('params','var')==0, error('Missing params struct'); end
if ~isfield(params,'lambda'), error('Missing lambda parameter'); end
assert(length(params.lambda)==1 || length(params.lambda) == nfreq, ...
       'lambda must be a scalar or a vector of length nfreq');

nlabs = parpool_size();
params.wri = true;
if isfield(params,'srcfreqmask'),srcfreqmask = params.srcfreqmask; else srcfreqmask = []; end
if nlabs==0
    [f,g,h] = PDEfunc(PDEopts.OBJ,m,Q,[],D,model,params,srcfreqmask);
else
    [f,g,h] = PDEfunc_dist(PDEopts.OBJ,m,Q,[],D,model,params,srcfreqmask);
end


end


