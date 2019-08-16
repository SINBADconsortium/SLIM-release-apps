function y = project(x,dimTree)
% PROJECT - Projects a general set of hierarchical tucker parameters x onto
% the set of orthonormalized parameters.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   y = project(x,dimTree)
%
% Input:
%   x       - vectorized, unorthogonalized parameters
%   dimTree - dimension tree object corresponding to x
% 
% Output:
%   y       - vectorized, orthogonalized parameters
    if length(size(x)) > 1 && min(size(x)) > 1 || length(x) ~= dimTree.sizeOfVec()
        error('invalid size');
    end
    [U,B] = dimTree.fromVec(x);
    T = dimTree.emptyTree(U,B);
    itr = dimTree.iterator('leaves',T);
    while itr.advance()
        U = itr.getValue('U');
        [Q,R] = qr_std(U);
        itr.setValue('U',Q);
        itr.setValue('R',R);
    end
    itr = dimTree.iterator('interior_up',itr.T);
    while itr.advance()
       Rleft = itr.getLeftValue('R');
       Rright = itr.getRightValue('R');
       B = itr.getValue('B');
       sz = size(B);
       if length(sz) == 2
           sz = [sz,1];
       end
       B = ttm(B,{Rleft,Rright},[1 2]);
       if ~itr.isRoot()
           P = matricize(B,[1 2],3);
           [Q,R] = qr_std(P);
           B = dematricize(Q,sz,[1 2],3);
           itr.setValue('R',R);
       end
       itr.setValue('B',B);
    end
    [U,B] = dimTree.extractParams(itr.T);
    y = dimTree.toVec(U,B);
           
end