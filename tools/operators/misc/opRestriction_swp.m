classdef opRestriction_swp < opSpot
    
    properties ( SetAccess = protected ) 
        idx;
    end
    
    methods 
        function op = opRestriction_swp(n,idx)
            op = op@opSpot('Restriction operator',length(idx),n);
            op.idx = idx; op.sweepflag = true;
        end
    end
    
    methods ( Access = protected )
        function y = multiply(op,x,mode)
            if mode == 1
                y = x(op.idx,:);
            else
                if iscodistributed(x) || isdistributed(x)
                    spmd,
                        xloc = getLocalPart(x); codist = getCodistributor(x);
                        yloc = zeros(size(op,2),size(xloc,2));
                        yloc(op.idx,:) = xloc;
                        codisty = codistributor1d(2,codist.Partition,[size(op,2),size(x,2)]);
                        y = codistributed.build(yloc,codisty,'noCommunication');
                    end
                else
                    y = zeros(size(op,2),size(x,2));
                    y(op.idx,:) = x;
                end
            end
        end
    end
end