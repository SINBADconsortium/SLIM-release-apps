function [f,g] = Jls(m,Q,D,model)
% least-squares misfit for frequency domain FWI
%
% use: 
%   [f,g] = Jls(m,Q,D,model)
%
% input:
%   m - gridded squared slowness [km^2/s^2], see also F.m
%   Q - gridded source functions, see also F.m
%   D - `observed' data consistent with output of F(m,Q,model)
%   model - struct with model parameters for F.m
%
% output
%   f - misfit
%   g - gradient
%

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

% model data
[Dt,Jt] = F(m,Q,model);

% misfit
f = .5*norm(Dt-D).^2;

% gradient, if needed
if nargout > 1
    g = Jt'*(Dt - D);
end
