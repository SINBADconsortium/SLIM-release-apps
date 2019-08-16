function [f,g] = JW(m,Q,D,W,model,params)
% misfit for simultaneous-source waveform inversion
%
% use:
%  [f,g] = JW(m,Q,D,W,model,params)
%
% input:
%   m - vector with gridded squared-slowness [km^2/s^2] (see also F)
%   Q - matrix containing gridded source functions (see also F)
%   D - observed data as distributed vector (consistent with output of F(m,Q,model))
%   W - sim. source matrix of size nsrc x nsimsrc
%   model - struct. with model parameters (see also F)
%   params.misfit  - function handle to penalty function of the form 
%                    [f,g] = misfit(x) where f is a scalar, g is the
%                    gradient of the misfit w.r.t. x 
%                    default: @(x)twonorms(x).
%
% output:
%   f - misfit value
%   g - gradient of misfit w.r.t. m

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if nargin < 6
    params = [];
end

% dimensions of data
nfreq = length(model.freq);
nrec  = length(model.xrec)*length(model.zrec);
nsrc  = size(Q,2);

% default options
misfit = @(x)twonorms(x);

% check options
if isfield(params,'misfit')
    misfit = params.misfit;
end

% sim. source 
WW  = oppBlockDiag(nfreq,opKron(W',opDirac(nrec)));
D   = WW*D;
Q   = Q*W;

% model data
[Dt,Jt] = F(m,Q,model);

% evaluate misfit
[f,df]  = misfit(Dt - D);

% gradient
if nargout > 1
    g  = Jt'*df;
end
