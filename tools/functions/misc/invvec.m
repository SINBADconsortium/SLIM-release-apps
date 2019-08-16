function x = invvec(x, sizes)
% Reverse of vec. Works like reshape but also for distributed arrays, given
% that it is distributed over the last dimension.
%
% use: 
%   function y = invvec(x, sizes)
% 
% input: 
%   x     - vector
%   sizes - required size
% 
% output:
%   y     - array of size [sizes].

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
    
    if not(isvector(x))
        error('distunvec: Input x needs to be a column vector')
    end
    
    assert(length(x) == prod(sizes), 'Length of input must equal product of sizes')
    nd = length(sizes);
    if nd > 1
        if isdistributed(x)
            x = pSPOT.utils.distVec2distArray(x,sizes);
        else
            x = reshape(x,sizes);
        end
    end
end
            