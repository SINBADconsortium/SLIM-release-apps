function H = opHalfInt(model)
%opHalfInt Returns operator that performs a trace by trace
% half-integration on a shot record of dimension ns x nrec.
%
%  H = opHalfInt(model)
%
%--------------------------------------------------------------------------
%
% INPUT:
%   model - structure containing the model parameters
%
% OUTPUT:
%   H  -  SPOT operator that performs the half integration on a
%         (vectorized) shot record
%
%

% Author: Philipp Witte (pwitte@eos.ubc.ca)
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: February 2016

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

    numSamples = length(model.NyqT);
    tmax = model.T;
    
    F1 = opDFTR(numSamples);
    dt = numSamples*1e3/tmax;
    df = dt/numSamples;
    freq = 0:df:dt/2;
    omega = freq*2*pi;
    partialDt = 1./sqrt(omega); partialDt(1) = 1;
    P1 = F1'*opDiag(partialDt)*F1;
    H = opKron(opDirac(length(model.xrec)),P1);

end