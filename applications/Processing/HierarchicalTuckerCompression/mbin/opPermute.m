classdef opPermute < opSpot
% OPPERMUTE - Permutation operator working on vectorized tensors. Working on
% distributed arrays requires the pSPOT toolbox.
%
% Curt Da Silva
% curtd@math.ubc.ca
%
% Usage:
%   P = opPermute(dims,permutation)
%
% Input:
%   dims            - dimensions of the input vectorized arrays
%   permutation     - permutation of 1:length(dims)
%
% Output:
%   P               - SPOT permutation operator
    
    properties(SetAccess = protected)
        dims, permutation;
    end
    
    methods 
        % dims - original dimensions of the vector
        function op = opPermute(dims, permutation)
            op = op@opSpot('Permutation',prod(dims),prod(dims));
            op.dims = dims;
            op.permutation = permutation;
            
        end
        function out = test(op)
            x = randn(prod(op.dims),1);
            y = randn(prod(op.dims),1);
            
            z = op.multiply(x,1);
            t = y' * z;
            
            z = op.multiply(y,-1);
            s = x' * z;
            
            e1 = abs(s - t)/abs(max(s,t));
                                         
            if~(e1<1e-10); fprintf(2,'opPermute: adjoint test failed, error = %g\n',e1); end
           
            z = op.multiply(x,1);
            y = op.multiply(z,-1);
            e2 = norm(x - y);
            if~(e2<1e-10); fprintf(2,'opPermute: not orthogonal,error = %f\n',e2);end                        
            
            out = (e1<1e-10) & (e2<1e-10);
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