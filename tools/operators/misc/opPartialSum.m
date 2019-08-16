classdef opPartialSum < opSpot
% opPartialSum - Partially sum a vector along a given set of indices
% 
% Curt Da Silva, 2016
%
% Usage:
%   op = opPartialSum(idx,sum_along_dim);
%
% Input:
%   idx           - [n x ndims] array of indices
%   sum_along_dim - dimension along which to sum
%
% For the opPartialSum operator P, the entries of the vector
%  x corresponds to multi-dimensional indices idx
% 
% y = P*x is such that
%  y(i) = sum(x(j) : such that j = i, for i \in idx(:,keep_dim) )
%       
    
    properties (SetAccess = protected)
        idx,keep_dim;
    end
    
    methods
        function op = opPartialSum(idx,keep_dim)
            m = length(unique(idx(:,keep_dim)));
            op = op@opSpot('opPartialSum',m,size(idx,1));
            op.idx = idx;
            op.keep_dim = keep_dim;
            op.sweepflag = true;
        end
    end
    
    methods (Access = protected)
        function y = multiply(op,x,mode)
            if mode==1                           
                [~,~,I] = unique(op.idx(:,op.keep_dim),'rows');
                p = size(x,2);
                m = size(op.idx,1);
                % Account for a number of right hand sides
                I = vec(repmat(I,1,p));
                rhs_idx = vec(repmat(1:p,m,1));
                I = sub2ind([size(op,1),p],I,rhs_idx);
                y = reshape(accumarray(I,vec(x)),size(op,1),p);
            else
                error('Not implemented yet');
            end
        end
    end
    
end