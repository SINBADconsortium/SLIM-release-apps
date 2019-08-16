function z = cell_skeleton(x)
% CELL_SKELETON - Recursively copies the cell structure of the input cell array
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca 
%
% Usage:
%   z = cell_skeleton(x);
%
% Input:
%   x - cell array, whose elements are cell arrays (up to a
%   particular depth)
%  
% Output:
%   z - cell array with the same cell structure as x
    
    z = {};
    if(~iscell(x))
        return;
    end

    for i=1:length(x)
        z{end+1} = cell_skeleton(x{i});
    end
end