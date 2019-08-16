function [f,g] = gramian_regularizer(dimTree, x, lambda)
% GRAMIAN_REGULARIZER - Gramian-based regularizer for objective
% functions defined on the HT manifold, i.e.
%
%  \lambda \sum_{t \in T} \|X^(t)\|F^2 + \|(X^{(t)})^{\dagger}\|F^2
%
% See 'Optimization on the Hierarchical Tucker manifold -
% applications to tensor completion', C. Da Silva and F. Herrmann,
% 2013, for more details.
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   [f,g] = gramian_regularizer(dimTree,x,lambda);
%
% Input:
%   dimTree - dimension tree object
%   x       - orthogonalized HT parameters
%   lambda  - weighting
% 
% Output:
%   f       - objective value above
%   g       - (optional) Riemannian gradient
    
G = dimTree.gramian(x);
T = dimTree.emptyTree();
T = dimTree.copyCell(T,G,'G');

f = 0;            
treeitr = dimTree.iterator('down',T);
epsilon = 1e-12;
while treeitr.advance()
    if ~treeitr.isRoot()
        [U,S,~] = svd(treeitr.getValue('G'));                    
        f = f + lambda * (sum(diag(S)) + sum(diag(S+epsilon).^(-1) ));
        if nargout == 2                        
            dG = lambda * ( U * diag(1 - diag(S+epsilon).^(-2)) * U');
            treeitr.setValue('dG',dG);
        end
    else
        treeitr.setValue('dG',0);
    end
        
end            

if nargout == 2
    dG = dimTree.extractParam(treeitr.T,'dG');
    dG = dimTree.toVecGram(dG);
    gramJ = opGramianJ(dimTree,x,G);
    g = gramJ' * dG;
end


end
