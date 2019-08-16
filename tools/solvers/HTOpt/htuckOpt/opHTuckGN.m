classdef opHTuckGN < opSpot
% OPHTUCKGN - Gauss-Newton Hessian 
% Equivalent to J' * J , where J is the HTucker Jacobian in opHTuckJ.m,
% but much much faster (all intermediate vectors are of the size of the
% tangent space, unlike naively applying J' * J). Can apply
% forward multiplication and inverse multiplication efficiently.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
%   Usage: 
%     H = opHTuckGN(dimTree,x)
% 
%   Input:
%     dimTree   - dimension tree object from which to take derivatives
%     x         - vectorized parameter matrices corresponding
%                 to the dimension tree (must be orthogonalized)
%   Output:
%     H         - SPOT operator for the Gauss-Newton Hessian
    properties(SetAccess = protected)
        dimTree,x, G;
    end
    
    methods        
        function op = opHTuckGN(dimTree,x)
            n = prod(dimTree.dims);
            op = op@opSpot('HT Gauss-Newton Hessian',length(dimTree.randn()), ...
                           length(dimTree.randn()));
            op.dimTree = dimTree;                       
            op.x = x;
            op.G = dimTree.gramian(op.x);
        end
        function y = mldivide(op, x)
            dimTree = op.dimTree;
            [dU,dB] = dimTree.fromVec(x);                                
            G = op.G;          
            T = dimTree.emptyTree(dU,dB,'dU','dB');
            T = dimTree.copyCell(T, G,'G');
            itr = dimTree.iterator('down',T);
            while itr.advance()
                if ~itr.isRoot()
                    G = itr.getValue('G');
                    if itr.isLeaf()
                        itr.setValue('dU',itr.getValue('dU')*opInverse(conj(G)));
                    else
                        itr.setValue('dB',ttm(itr.getValue('dB'),opInverse(G),3));
                    end
                end
            end
            [dU,dB] = dimTree.extractParams(itr.T,'dU','dB');
            y = dimTree.toVec(dU,dB);
        end
    end
    
    methods(Access = protected)
        
        function y = multiply(op,x,mode)
            dimTree = op.dimTree;
            [dU,dB] = dimTree.fromVec(x);                                
            G = op.G;          
            T = dimTree.emptyTree(dU,dB,'dU','dB');
            T = dimTree.copyCell(T, G,'G');
            itr = dimTree.iterator('down',T);
            while itr.advance()
                if ~itr.isRoot()
                    G = itr.getValue('G');
                    if itr.isLeaf()
                        itr.setValue('dU',itr.getValue('dU')*conj(G));
                    else
                        itr.setValue('dB',ttm(itr.getValue('dB'),G,3));
                    end
                end
            end
            [dU,dB] = dimTree.extractParams(itr.T,'dU','dB');
            y = dimTree.toVec(dU,dB);
            
        end

    end
end