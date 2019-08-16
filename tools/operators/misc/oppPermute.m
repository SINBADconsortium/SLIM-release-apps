classdef oppPermute < oppSpot
% OPPPERMUTE - Permutation operator working on vectorized,
% distributed tensors. 
%
% Curt Da Silva, 2014
% curtd@math.ubc.ca
%
% Usage:
%   P = oppPermute(dims,permutation)
%
% Input:
%   dims            - dimensions of the input vectorized arrays
%   permutation     - permutation of 1:length(dims)
%
% Output:
%   P               - pSPOT permutation operator
    
    properties(SetAccess = protected)
        dims, permutation;
    end
    
    methods 
        % dims - original dimensions of the vector
        function op = oppPermute(dims, permutation)
            op = op@oppSpot('Permutation',prod(dims),prod(dims));
            op.dims = dims;
            op.permutation = permutation;
            
        end
    end
    
    methods(Access = protected)
        function y = multiply(op,x,mode)            
            if(mode == 1)
                p = op.permutation;
                dims = op.dims;
            else
                dims = op.dims(op.permutation);
                p(op.permutation) = 1:length(op.permutation);
            end            
            
            if(isdistributed(x))
                y = pSPOT.utils.distVec2distArray(x,dims);
                y = pSPOT.utils.distVectorize(pSPOT.utils.DistPermute(y,p));
            else
                y = vec(permute(reshape(x,dims), p));
            end
         
        end
    end
    
end