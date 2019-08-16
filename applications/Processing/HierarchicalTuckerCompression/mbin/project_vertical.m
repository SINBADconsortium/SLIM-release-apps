function z = project_vertical(x,dx,dimTree)
% PROJECT_VERTICAL - Extracts the component of the tangent vector in the vertical space
% i.e. in the nullspace of opHTuckJ    
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   z = project_vertical(x,dx,dimTree)
%
% Input:
%   x   - (orthogonalized) parameters
%   dx  - tangent vector. length(dx) == length(x)
%
% Output:
%   z   - vertical component of dx
    
[U,B] = dimTree.fromVec(x);
[dU,dB] = dimTree.fromVec(dx);

T = dimTree.emptyTree(U,B);
T = dimTree.copyParams(T,dU,dB,'dU','dB');
itr = dimTree.iterator('up',T);
skew = @(A) (A + A')/2;
while itr.advance()
   if itr.isLeaf()
       U = itr.getValue('U'); dU = itr.getValue('dU');
       D = skew(U' * dU);
       itr.setValue('D',D);
       dU = U * D;
       itr.setValue('dU',dU);
   else
       B = itr.getValue('B'); dB = itr.getValue('dB');
       Z = ttm(B,itr.getLeftValue('D'),1) + ttm(B,itr.getRightValue('D'),2);
       if itr.isRoot()
          dB = -Z; 
       else          
          D = skew(ttt(B,dB + Z,[1 2],[1 2]));          
          dB = ttm(B,D',3) - Z;
          itr.setValue('D',D);
       end
       itr.setValue('dB',dB);
   end
end

[dU,dB] = dimTree.extractParams(itr.T,'dU','dB');
z = dimTree.toVec(dU,dB);

end