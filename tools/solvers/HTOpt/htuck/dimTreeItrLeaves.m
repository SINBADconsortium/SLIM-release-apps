classdef dimTreeItrLeaves < handle & dimTreeItr
% DIMTREEITERATOR - iterator class for dimensionTree objects and accompanying 
% cell arrays which have the same structure. These objects should not be
% initialized directly, but should be generated from an existing dimTree
% object.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
%
% Usage:
%   itr = dimTree.iterator(mode,T);
%   while ~itr.isFinished()
%     %CODE
%     itr.advance();
%   end
%   results = itr.T;
%
% See dimTree.iterator for documentation and proper initialization
%
% Available methods:
%
%%%%% Navigation
%   isFinished()                   - returns true if the iterator has finished iterating, false otherwise
%   advance()                      - advances the iterator to its next position
%   isEndOfRow()                   - returns true if the iterator is at the end of the current level, false otherwise
%   isRoot()                       - returns true if the iterator is at the root node, false otherwise
%   isLeaf()                       - returns true if the iterator is at a leaf node, false otherwise
%
%%%%% Data storage/retrieval
%   getValue(var_name)             - returns the value of the variable named var_name associated to the current node
%   setValue(var_name,var_val)     - sets the value var_val of the variable named var_name to the current node
%   getLeftValue(var_name),        - same as getValue, but for the left and right children, respectively, if 
%   getRightValue(var_name)          not at a leaf node
%   setLeftValue(var_name,var_val) - same as setValue, but for the left and right children, respectively, if 
%   setRightValue(var_name,var_val)  not at a leaf node
%   getLabels()                    - returns the names of all the variables associated to the current node
%   getDims()                      - get the dimensions associated to the current node
%   getRank()                      - get the rank associated to the current node
       
    
    methods
        function itr = dimTreeItrLeaves(dimTree,T)
        % Construct a dimTreeIterator object. Should not be called
        % directly, use dimTree.iterator(mode,T) instead.
        %
        % Usage:
        %    itr = leavesIterator(dimTree);
        %
        % Input:
        %    dimTree - dimensionTree object specifying the tree structure
        %              to iterate over
        %    mode    - up: iterates from leaves to root        
        %              down: iterates from root to leaves
        %              interior_up: iterates from base of tree (ignoring
        %              leaves) to root
        %              interior_down: iterates from root to base of tree,
        %              ignoring leaves
        %              leaves: iterates over leaves only
        % 
        %   T        - cell array of structs conforming to dimTree. Used to
        %              store/retrieve values associated to each node.
            itr = itr@dimTreeItr(dimTree,T);
        end       
        function z = getValue(itr,name)
        % Gets the value associated to the current node.
        %
        % Usage: 
        %   z = itr.getValue(name)
        %
        % Input:
        %   name - variable name associated to the current node
        %
        % Output:
        %   z    - variable with the provided input name
            [level,idx] = itr.leafIdx2d();
            z = getValue@dimTreeItr(itr,name,level,idx);
        end 
        
        function setValue(itr,name,z)
        % Sets the value associated to the current node.
        %
        % Usage: 
        %   itr.setValue(name,z)
        %
        % Input:
        %   name - variable name associated to the current node
        %   z - variable to store        
            [level,idx] = itr.leafIdx2d();
            setValue@dimTreeItr(itr,name,z,level,idx);
        end                
        
        function labels = getLabels(itr)
        % Returns a cell array of all of the variable names associated to
        % the current node.
        % 
        % Usage:
        %   labels = itr.getLabels();
        %  
        % Output:
        %   labels - cell array of variable names at the current node
            [level,idx] = itr.leafIdx2d();
            labels = fieldnames(itr.T{level}{idx});
        end
        function dims = getDims(itr)
        % Returns the dimensions associated to the current node of the dimension tree
        % 
        % Usage:
        %   dims = itr.getDims();
        % 
        % Output:
        %   dims - vector of dimensions associated to the current node
            dims = itr.idx;
        end        
        
        function rank = getRank(itr)
        % Returns the rank associated to the current node of the dimension tree.
        %
        % Usage:
        %   rank = itr.getRank();
        %            
            [level,idx] = itr.leafIdx2d();
            rank = itr.dimTree.rank(level,idx);             
        end
        function z = advance(itr)
            if 0 <= itr.idx && itr.idx < itr.dimTree.d                
                itr.idx = itr.idx + 1;
                z = true;
            else 
                itr.idx = -1;
                z = false;
            end
        end        
    end       
    
    methods (Access = protected)
        function [level,idx] = leafIdx2d(itr)   
            level = itr.dimTree.leaf_indices(1,itr.idx);
            idx = itr.dimTree.leaf_indices(2,itr.idx); 
        end                       
    end
    
end

