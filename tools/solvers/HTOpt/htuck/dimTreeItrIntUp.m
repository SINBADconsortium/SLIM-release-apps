classdef dimTreeItrIntUp < handle & dimTreeItr
% DIMTREEUPITERATOR - iterator class for dimensionTree objects and accompanying 
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
        function itr = dimTreeItrIntUp(dimTree,T)
            itr = itr@dimTreeItr(dimTree,T);            
        end
        function z = advance(itr)
            if ~itr.isFinished()
                if itr.idx == 0
                    itr.depth = itr.dimTree.depth-1;
                    itr.idx = 1;
                elseif itr.idx == length( itr.dimTree.nodes{itr.depth} )
                    if itr.depth == 1
                        itr.idx = -1;
                        z = false; return;
                    else
                        itr.depth = itr.depth - 1; itr.idx = 1;
                    end
                else
                    itr.idx = itr.idx + 1;                  
                end 
                %Advance to next non-leaf node
                while itr.isLeaf()
                    if itr.idx == length( itr.dimTree.nodes{itr.depth} )
                        if itr.depth == 1
                            itr.idx = -1; z = false; return;
                        else
                            itr.depth = itr.depth - 1; itr.idx = 1;
                        end
                    else
                        itr.idx = itr.idx + 1;
                    end
                end
                z = true;
            else
                z = false;
            end
        end                                       
    end                                          
end

