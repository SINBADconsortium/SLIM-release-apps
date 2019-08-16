function y = project_horizontal(x, tvec, dimTree, zeroOutRoot)
% PROJECT_HORIZONTAL - Project a tangent vector onto the horizontal space.
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%    y = project_horizontal(x,tvec,dimTree,zeroOutRoot)
%
% Input:
%   x           - current point on the manifold (assumed to be orthogonalized)
%   tvec        - tangent vector
%   zeroOutRoot - true: sets the tangent vector at the root to be zero (default: false)
%
% Output:
%   y           - horizontal vector
    if(exist('zeroOutRoot','var')==0)
        zeroOutRoot = false;
    end
    
    [U,B] = dimTree.fromVec(x);
    T = dimTree.emptyTree(U,B);
    
    [dU,dB] = dimTree.fromVec(tvec);
    T = dimTree.copyParams(T,dU,dB,'dU','dB');
    
    itr = dimTree.iterator('down',T);
    while itr.advance()
       if ~itr.isRoot()
           if itr.isLeaf()
               U = itr.getValue('U');
               dU = itr.getValue('dU');
               dU = opProjection(U,true,true) * dU;
               itr.setValue('dU',dU);
           else
               B = itr.getValue('B'); dB = itr.getValue('dB');
               dB = opProjection(matricize(B,[1 2]),true,true) * matricize(dB,[1 2]);
               itr.setValue('dB',dematricize(dB,size(B),[1 2]));
           end
       elseif zeroOutRoot
           itr.setValue('dB',zeros(size(itr.getValue('dB'))));
       end
    end
    [dU,dB] = dimTree.extractParams(itr.T,'dU','dB');    
    y = dimTree.toVec(dU,dB);               

end
