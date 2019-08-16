classdef opHTuckJ2 < opSpot
% OPHTUCKJ2 - Jacobian of the function the full-tensor expansion operator 
% (i.e. dimTree.full(x)) - serial version. This Jacobian uses orthogonal
% projection to ensure that the input/output vectors are in the horizontal space 
% at the current point. Avoids computing P_{U_t}^{\perp} when t is
% not a leaf node, which is identical to opHTuckJ mathematically
% but cheaper to compute.
%
% See oppHTuckJ for a parallel version of this code.
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Use:
%   J = opHTuckJ2(dimTree,x,{fullTree})
%   
% Input:
%   dimTree          - dimension tree object from which to take derivatives
%   x                - vectorized parameter matrices corresponding
%                      to the dimension tree
%   fullTree         - if not empty, the result of dimTree.fullTree(), 
%                      if previously computed elsewhere (otherwise,
%                      it is computed in this method) 

    properties(SetAccess = protected)
        dimTree,fullTree,x;
    end
    methods        
        function op = opHTuckJ2(dimTree,x, fullTree)
            n = prod(dimTree.dims);
            op = op@opSpot('Hierarchical Tucker Jacobian',n, ...
                           length(dimTree.randn()));
            op.dimTree = dimTree;            
            if(exist('fullTree','var') && ~isempty(fullTree))
                op.fullTree = fullTree;
            else                
                [U,B] = dimTree.fromVec(x);
                op.fullTree = dimTree.fullTree(U,B);
            end            
            op.x = x;
        end
    end
    methods(Access = protected)
        function y = multiply(op,x,mode)                        
            dimTree = op.dimTree;
            fullTree = op.fullTree;                       
            
            if(mode == 1)
                [dU,dB] = dimTree.fromVec(x);
                T = dimTree.emptyTree(dU,dB,'dU','dB');
                T = dimTree.copyTree(T, fullTree, @(x) strcmp(x,'U') || strcmp(x,'B'));
                itr = dimTree.iterator('leaves',T);
                while itr.advance()
                   U = itr.getValue('U');
                   dU = itr.getValue('dU');                   
                   itr.setValue('dU',dU);
                end
                
                itr = dimTree.iterator('interior_up',itr.T);
                while itr.advance()
                   UL = itr.getLeftValue('U'); UR = itr.getRightValue('U'); B = itr.getValue('B');
                   dUL = itr.getLeftValue('dU'); dUR = itr.getRightValue('dU'); dB = itr.getValue('dB');                   
                   dU = ttm(B,{dUL,UR},[1 2]) + ttm(B,{UL,dUR},[1 2]) + ttm(dB,{UL,UR},[1 2]);                     
                   dU = matricize(dU,[1 2]);                   
                   itr.setValue('dU',dU);
                end
                y = itr.T{1}{1}.dU;                                
            else
                T = dimTree.emptyTree();
                T = dimTree.copyTree(T,fullTree, @(x) strcmp(x,'U') || strcmp(x,'B'));
                T{1}{1}.dU = x;
                
                itr = dimTree.iterator('interior_down', T);
                while itr.advance()
                   UL = itr.getLeftValue('U'); UR = itr.getRightValue('U'); B = itr.getValue('B');
                   dU = itr.getValue('dU');
                   
                   nl = size(UL,1); nr = size(UR,1); k = size(B,3);
                   if itr.isRoot()
                       if(issparse(dU))
                           dU = full(dU);
                       end
                       dU = reshape(dU,nl,nr);
                       dUL1 = dU * conj(UR);              
                       dUL = ( dUL1 * B' );
                       dUR = ( ((UL * B)' * dU ) .' ) ;
                       dB = UL' * dUL1;    
                       clear dUL1 x;
                   else
                       dU = dematricize(dU,[nl, nr, k],[1 2]);
                       dUL = (ttt(B,ttm(dU,UR',2),[2 3],[2 3]) .');
                       dUR = (ttt(B,ttm(dU,UL',1),[1 3],[1 3]) .');
                       dB = ttm(dU,{UL',UR'},[1 2]);
                      
                   end
                   itr.setValue('dB',dB); itr.setLeftValue('dU',dUL); itr.setRightValue('dU',dUR);                       
                   clear dU dUL dUR dB UL UR B;
                end               
                
                [dU,dB] = dimTree.extractParams(itr.T,'dU','dB');

                y = dimTree.toVec(dU,dB);
                y = project_horizontal(op.x,y,dimTree);
                
            end
        end
    end
end
