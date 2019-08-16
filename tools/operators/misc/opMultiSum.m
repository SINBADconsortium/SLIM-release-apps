classdef opMultiSum < opSpot
    
    properties (SetAccess = protected)
        ops;
    end
    methods
        function op = opMultiSum(operators)
            assert(isa(operators,'cell') && length(find(cellfun(@(x) isa(x,'opSpot'),operators)))==length(operators),'operators must be a cell array of spot operators');
            assert(max(cellfun(@(x) size(x,1),operators))==size(operators{1},1) && max(cellfun(@(x) size(x,2),operators))==size(operators{1},2),'All SPOT operators must have the same size');
            op = op@opSpot('Sum of SPOT operators',size(operators{1},1),size(operators{1},2));
            op.ops = operators;
        end        
    end
    methods (Access = protected)        
        function y = multiply(op,x,mode)
            if mode==1
                y = zeros(size(op,1),1);
                for i=1:length(op.ops)
                    y = y + op.ops{i}*x;
                end
            else
                y = zeros(size(op,2),1);
                for i=1:length(op.ops)
                    y = y + op.ops{i}'*x;
                end
            end
        end
    end
end