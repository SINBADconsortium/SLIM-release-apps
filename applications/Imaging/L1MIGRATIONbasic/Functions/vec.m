function output = vec(input)
% Vectorizes input. Use invvec to undo this operation.
%
% use:
%   output = vec(input);
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

if isdistributed(input)
    output = pSPOT.utils.distVectorize(input);
else
    output = input(:);
end