function x = treeCell2Vec(C)
% TREECELL2VEC - Vectorizes a cell array (arranged as a tree) of
% matrices/tensors, from the root to the leaves.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
% 
% Usage:
%   x - treeCell2Vec(C)
% 
% Input:
%   C - cell array of matrices/tensors
% 
% Output:
%   x - vectorized cell array
   total = 0;
   for i=1:length(C)
       for j=1:length(C{i})
           total = total + numel(C{i}{j});
       end
   end
   x = zeros(total,1);
   total = 1;
   for i=1:length(C)
       for j=1:length(C{i})
           x(total:total+numel(C{i}{j})-1) = vec(C{i}{j});
           total = total + numel(C{i}{j});
       end
   end
end