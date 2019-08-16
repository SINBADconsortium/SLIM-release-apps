function s = getoption(options,field,default)
% parse option from struct
%
% use
%   s = getoption(options,field,default)
%
% input
%   options - struct with option fields
%   field   - specify which field to parse
%   default - default value in case field does not exist
%
% output
%   value
%
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

s = default;
if isfield(options,field)
    s = getfield(options,field);
end