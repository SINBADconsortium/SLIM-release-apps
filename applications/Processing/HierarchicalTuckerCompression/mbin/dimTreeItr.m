classdef (Abstract,HandleCompatible) dimTreeItr
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
    
    properties
        dimTree,T,depth,idx;
    end
    
    methods (Abstract)
        advance;         
    end
    
    methods
        function itr = dimTreeItr(dimTree,T)
            itr.dimTree = dimTree;
            itr.T = T;
            itr.depth = 0;
            itr.idx = 0;
        end
        function z = isFinished(itr)
            z = itr.idx == -1;
        end
        
        function z = isRoot(itr)
            z = (itr.depth == 1 && itr.idx == 1);
        end
        function z = getValue(itr,name, level,idx)
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
            if nargin == 2
                z = itr.getValueIdx(name);                
            else
                z = itr.getValueIdx(name,level,idx);
            end
        end 
        function setValue(itr,name,z,level,idx)
        % Sets the value associated to the current node.
        %
        % Usage: 
        %   itr.setValue(name,z)
        %
        % Input:
        %   name - variable name associated to the current node
        %   z - variable to store        
            if nargin == 3
                itr.setValueIdx(name,z);
            else
                itr.setValueIdx(name,z,level,idx);
            end
        end
        
        function z = isLeaf(itr)
        % Returns true if the iterator is at a leaf node, false otherwise.
            z = length(itr.dimTree.nodes{itr.depth}{itr.idx}.dims)==1;
        end
        function z = getLeftValue(itr,name)
        % Gets the value associated to the left child of the current node.
        %
        % Usage: 
        %   z = itr.getLeftValue(name)
        %
        % Input:
        %   name - variable name associated to the left child of the current node
        %
        % Output:
        %   z    - variable with the provided input name    
            if ~itr.isLeaf()
                left_idx = 2*(itr.idx-1)+1;
                z = itr.getValueIdx(name,itr.depth+1,left_idx);
            end
        end
        function z = getRightValue(itr,name)
        % Gets the value associated to the right child of the current node.
        %
        % Usage: 
        %   z = itr.getLeftValue(name)
        %
        % Input:
        %   name - variable name associated to the right child of the current node
        %
        % Output:
        %   z    - variable with the provided input name        
            if ~itr.isLeaf()                
                right_idx = 2*(itr.idx-1)+2;
                z = itr.getValueIdx(name,itr.depth+1,right_idx);
            end
        end
        function setLeftValue(itr,name,z)
        % Sets the value associated to the left child of the current node.
        %
        % Usage: 
        %   itr.setLeftValue(name,z)
        %
        % Input:
        %   name - variable name associated to the left child of the current node
        %   z - variable to store      
            if itr.isLeaf()
                error('At a leaf node, no left child exists');
            else
                left_idx = 2*(itr.idx-1)+1;
                itr.setValueIdx(name,z,itr.depth+1,left_idx);
            end                     
        end
        function setRightValue(itr,name,z)
        % Sets the value associated to the right child of the current node.
        %
        % Usage: 
        %   itr.setRightValue(name,z)
        %
        % Input:
        %   name - variable name associated to the right child of the current node
        %   z - variable to store 
            if ~itr.isLeaf()
                right_idx = 2*(itr.idx-1)+2;
                itr.setValueIdx(name,z,itr.depth+1,right_idx);
            end              
        end
         function labels = getLabels(itr)
        % Returns a cell array of all of the variable names associated to
        % the current node.
        % 
        % Usage:
        %   labels = itr.getLabels();
        %  
        % Output:
        %   labels - cell array of variable names at the current
        %   node
            labels = fieldnames(itr.T{itr.depth}{itr.idx});
        end
        function dims = getDims(itr)
        % Returns the dimensions associated to the current node of the dimension tree
        % 
        % Usage:
        %   dims = itr.getDims();
        % 
        % Output:
        %   dims - vector of dimensions associated to the current node
            dims = itr.dimTree.nodes{itr.depth}{itr.idx}.dims;
        end
        function dims = getLeftDims(itr)
            if ~itr.isLeaf()
                left_idx = 2*(itr.idx-1)+1;
                dims = itr.dimTree.nodes{itr.depth+1}{left_idx}.dims;
            end
        end
        function dims = getRightDims(itr)
            if ~itr.isLeaf()
                right_idx = 2*(itr.idx-1)+2;
                dims = itr.dimTree.nodes{itr.depth+1}{right_idx}.dims;
            end
        end
        function rank = getRank(itr)
        % Returns the rank associated to the current node of the dimension tree.
        %
        % Usage:
        %   rank = itr.getRank();
        %
            rank = itr.dimTree.nodes{itr.depth}{itr.idx}.rank;
        end
        function rank = getLeftRank(itr)
            if ~itr.isLeaf()
                left_idx = 2*(itr.idx-1)+1;
                rank = itr.dimTree.nodes{itr.depth+1}{left_idx}.rank;
            end
        end
        function rank = getRightRank(itr)
            if ~itr.isLeaf()
                right_idx = 2*(itr.idx-1)+2;
                rank = itr.dimTree.nodes{itr.depth+1}{right_idx}.rank;
            end
        end
    end                                    
    
    methods (Access = protected)       
        function z = getValueIdx(itr,name,depth,idx)
            if itr.isFinished()
                error('Tried to get value past the iterator end');
            end
            if exist('name','var')==0 || isempty(name)
                error('Need to specify field');
            end
            if nargin == 2
                depth = itr.depth; idx = itr.idx;
            end
            z = itr.T{depth}{idx}.(name);
        end
      
        function setValueIdx(itr,name,z,depth,idx)
            if itr.isFinished()
                error('Tried to set value past the iterator end');
            end                     
            if nargin == 3
                depth = itr.depth; idx = itr.idx;
            end                                
            itr.T{depth}{idx}.(name) = z;    
        end
    end
    
end

