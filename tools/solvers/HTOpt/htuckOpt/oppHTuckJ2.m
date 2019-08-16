classdef oppHTuckJ2 < oppSpot
% OPPHTUCKJ - Jacobian of the function the full-tensor expansion operator 
% (i.e. dimTree.full(x)) - parallel version. This Jacobian uses orthogonal
% projection to ensure that the input/output vectors are in the horizontal space 
% at the current point.
%
% This code acts on parallel vectors in the (full) ambient space,
% it is serial in the non-root nodes of the dimension tree.
%
% See opHTuckJ for a serial version of this code
%
% Curt Da Silva
% HTOpt v1.0
% curtd@math.ubc.ca
%
% Use:
%   J = oppHTuckJ(dimTree,x,{fullTree})
%   
% Input:
%   dimTree          - dimension tree object from which to take derivatives
%   x                - vectorized parameter matrices corresponding 
%                    to the dimension tree
%   fullTree         - if not empty, the result of dimTree.fullTree(),
%                    if previously computed elsewhere (otherwise, 
%                    it is computed in this method) (default: [])

    properties(SetAccess = protected)
        dimTree,fullTree,x;
    end
    methods
        function op = oppHTuckJ2(dimTree,x, fullTree)
            n = prod(dimTree.dims);
            op = op@oppSpot('Parallel Hierarchical Tucker Jacobian',n, ...
                           length(dimTree.randn()));
            op.dimTree = dimTree;            
            if(exist('fullTree','var') && ~isempty(fullTree))
                op.fullTree = fullTree;
            else                
                [U,B] = dimTree.fromVec(x);
                op.fullTree = dimTree.fullTree(U,B,false);
            end
            
            op.x = x;
        end
    end
    methods(Access = protected)
        function y = multiply(op,x,mode)
            derivativeTree = cell_skeleton(op.fullTree);
            %if forward mode, derivatives get sent up the tree
            %if reverse mode, derivatives get sent down the tree
            dimTree = op.dimTree;
            fullTree = op.fullTree;
            if(mode == 1)
                error('Not implemented for parallel');
            else                
                T = dimTree.emptyTree();
                T = dimTree.copyTree(T,fullTree, @(x) strcmp(x,'U') || strcmp(x,'B'));                
                
                itr = dimTree.iterator('interior_down', T);
                while itr.advance()
                   UL = itr.getLeftValue('U'); UR = itr.getRightValue('U'); B = itr.getValue('B');                   
                   
                   nl = size(UL,1); [nr,kr] = size(UR); k = size(B,3);
                   gather = 1;
                   if itr.isRoot()
                       dUL1 = reshape( oppKron2Lo(opMatrix(UR'),opDirac(nl),gather) * x,[nl,kr]);
                       dUL = ( dUL1 * B' );
                       dUR = oppKron2Lo(opDirac(nr), opMatrix((UL * B)'),gather) * x;
                       dUR = reshape(opPermute([kr,nr],[2 1]) * dUR,[nr,kr]);
                       dB = UL' * dUL1;                      
                   else
                       dU = itr.getValue('dU');
                       dU = dematricize(dU,[nl, nr, k],[1 2]);
                       dUL = ( ttt(B,ttm(dU,UR',2),[2 3],[2 3]) .');
                       dUR = ( ttt(B,ttm(dU,UL',1),[1 3],[1 3]) .');
                       dB = ttm(dU,{UL',UR'},[1 2]);
                   end
                   itr.setValue('dB',dB); itr.setLeftValue('dU',dUL); itr.setRightValue('dU',dUR);
                end
                [dU,dB] = dimTree.extractParams(itr.T,'dU','dB');
                y = dimTree.toVec(dU,dB);                                
                y = project_horizontal(op.x,y,dimTree);
            end
        end
    end
end