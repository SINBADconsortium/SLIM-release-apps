function A = vec2TreeCell(x, tree)
% VEC2TREECELL - Unvectorizes a vectorized cell array (assumed to
% represent a tree).
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%  A = vec2TreeCell(x,tree);
%
% Input:
%   x    - vectorized parameters
%   tree - cell array, whose structure A will take on, containing a
%          size vector (of integers) corresponding to the size of the
%          matrix/tensor at that node
%
% Output:
%   A    - tree cell array containing matrices/tensors
   A = cell_skeleton(tree);
   total = 1;
   for i=1:length(tree)
       for j=1:length(tree{i})
           if(~isempty(tree{i}{j}))
               dims = tree{i}{j};
               A{i}{j} = reshape(x(total:total+prod(dims)-1),dims);
               total = total+prod(dims);
           end
       end
   end
end