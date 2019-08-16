classdef dimensionTree < handle
% DIMENSION_TREE - Object representing a dimension tree for the hierarchical tucker format.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% 
% Usage:
%   dimTree = dimensionTree(dims, max_leaf_rank, max_interior_rank)
%   
% Input: 
%   dims              - size of each dimension 
%   max_leaf_rank     - rank of each leaf node (corresponding to
%                       the rank of the matricizations of each singleton dimension)
%   max_interior_rank - rank of each interior node (corresponding
%                       to the rank of the matricizations for each interior indices)
    properties 
       d, dims, root, nodes, depth, max_leaf_rank, max_interior_rank, leaf_indices
    end
    
    methods
        function dimTree = dimensionTree(dims, max_leaf_rank, max_interior_rank)
        % Constructs a balanced, complete HT dimension tree.
        % 
        % dimTree = dimensionTree(dims,max_leaf_rank, max_interior_rank);
        %
        % Input:
        %   dims - dimensions of the full tensor
        %   max_leaf_rank - rank associated to each leaf node
        %   max_interior_rank - rank associated to each internal node
        %
        % Output:
        %   dimTree - dimension tree object
            d = length(dims);
            depth = ceil(log2(d))+1;
            dimTree.root = struct;
            dimTree.root.size = dims;
            dimTree.root.dims = 1:d;
            dimTree.root.rank = 1;
            
            dimTree.dims = dims;
            dimTree.d = d;
            dimTree.depth = depth;
            dimTree.nodes = {{}};
            dimTree.max_leaf_rank = max_leaf_rank;
            dimTree.max_interior_rank = max_interior_rank;
            
            dimTree.nodes = cell(1,depth);
            dimTree.nodes{1} = {dimTree.root};
            leaf_indices = zeros(2,d);
            for i=2:depth
                if i < depth
                    nnodes = 2^(i-1);                    
                else
                    nnodes = d;
                    for j=1:length(dimTree.nodes{i-1})
                        if ~isempty(dimTree.nodes{i-1}{j}) && length(dimTree.nodes{i-1}{j}.dims)==1
                            nnodes = nnodes - 1;
                        end
                    end
                end                
                dimTree.nodes{i} = cell(1,nnodes);
                for j=1:2:nnodes                    
                    parent_idx = floor((j-1)/2)+1;
                    parent_dims = dimTree.nodes{i-1}{parent_idx}.dims;
                    p = length(parent_dims);
                    if p==1
                        continue;
                    end
                    if 2^(ceil(log2(p)))==p
                        mid = p/2;
                    else
                        mid =  ceil(p/2);
                        while ~( 2^(ceil(log2(mid))) == mid || 2^(ceil(log2(p-mid))) == p-mid )
                            mid = mid+1;
                        end
                    end
                    left_dims = parent_dims(1:mid);
                    right_dims = parent_dims(mid+1:end);
                    left_node = struct; right_node = struct;
                    left_node.dims = left_dims;
                    right_node.dims = right_dims;
                    if length(left_dims)==1
                        left_node.rank = max_leaf_rank;
                        leaf_indices(:,left_dims) = [i; j]; 
                    else
                        left_node.rank = max_interior_rank;
                    end
                    if length(right_dims)==1
                        right_node.rank = max_leaf_rank;
                        leaf_indices(:,right_dims) = [i; j+1];
                    else
                        right_node.rank = max_interior_rank;
                    end
                    dimTree.nodes{i}{j} = left_node;
                    dimTree.nodes{i}{j+1} = right_node;                    
                end
            end            
            dimTree.leaf_indices = leaf_indices;
        end
        
        function [U,B] = fromVec(dimTree, v)
        % Reshapes a vectorized set of parameter matrices into its
        % cell-array form
        %
        % Use:
        %   [U,B] = dimTree.fromVec(x)
        % 
        % Input: 
        %   x - vectorized HT parameters
        % 
        % Output:
        %   U - HT leaf bases
        %   B - HT interior nodes
            if length(size(v)) > 1 && (min(size(v)) > 1) || length(v) ~= dimTree.sizeOfVec()
                error('invalid size');
            end
            d = dimTree.d;
            dims = dimTree.dims;
            U = cell(d,1);
            total = 1;
            for i=1:d     
                idx = dimTree.leaf_indices(:,i);
                dim_size = dims(i);
                len = dims(i) * dimTree.nodes{idx(1)}{idx(2)}.rank;
                U{i} = reshape(v(total:total+len-1),dims(i),dimTree.nodes{idx(1)}{idx(2)}.rank);
                total = total + len;
            end            
            B = cell(dimTree.depth-1,1);            
            
            for i=1:dimTree.depth-1
                B{i} = cell(length(dimTree.nodes{i}),1);
                for j=1:length(B{i})
                    if length(dimTree.nodes{i}{j}.dims)==1
                        B{i}{j} = [];
                    else
                        left_idx = 2*(j-1)+1; right_idx = 2*(j-1)+2;
                        dims = [dimTree.nodes{i+1}{left_idx}.rank, ...
                                dimTree.nodes{i+1}{right_idx}.rank,...
                                dimTree.nodes{i}{j}.rank];
                        len = prod(dims);
                        B{i}{j} = reshape(v(total:total+len-1),dims);
                        total = total + len;
                    end                                       
                end                                
            end            
        end
        
        function T = emptyTree(dimTree, U,B, Uname, Bname)
        % Generates a cell array tree with the same structure as the
        % current dimension tree. 
        % 
        % Usage:
        %    T = dimTree.emptyTree({U},{B},{Uname},{Bname})
        %
        % Input:
        %    U     - cell array of variables associated to the leaves (optional)
        %    B     - cell array of variables associated to the interior nodes (optional)
        %    Uname - variable name for the leaf variables (default: 'U')
        %    Bname - variable name for the interior variables (default: 'B')
        %
        % Output:
        %    T     - cell array of structs corresponding to dimTree, with
        %            the specified variables at each node, if provided
           T = cell_skeleton(dimTree.nodes); 
           for i=1:length(T)
               for j=1:length(T{i})
                   T{i}{j} = struct;
               end
           end
           if exist('U','var') && exist('B','var')
               if exist('Uname','var')==0
                   Uname = 'U';
               end
               if exist('Bname','var')==0
                   Bname = 'B';
               end
               T = dimTree.copyParams(T,U,B, Uname, Bname);
           end
        end
        function T = copyParams(dimTree,T,U,B, Uname, Bname)
        % Write leaf/interior node parameters to the input cell-array.
        %
        % Usage:
        %   T = dimTree.copyParams(T,U,B,{Uname},{Bname});
        % 
        % Input:
        %   T     - cell array of structs corresponding to dimTree
        %   U     - cell array of variables associated to the leaves
        %   B     - cell array of variables associated to the interior nodes
        %   Uname - variable name associated to the leaves (default: 'U')
        %   Bname - variable name associated to the interior nodes (default: 'B')            
            if exist('Uname','var')==0
               Uname = 'U'; 
            end
            if exist('Bname','var')==0
               Bname = 'B'; 
            end
            itr = dimTree.iterator('down',T);
            while itr.advance()
                if itr.isLeaf()
                    leaf = dimTree.nodes{itr.depth}{itr.idx};
                    leafCtr = leaf.dims;
                    itr.setValue(Uname,U{leafCtr});
                else
                    itr.setValue(Bname,B{itr.depth}{itr.idx});
                end
            end
            
            T = itr.T;
        end
        function itr = iterator(dimTree,mode,T)
        % Generates an iterator corresponding to this dimensionTree.
        %
        % Usage:
        %    itr = dimTree.iterator(mode,T);
        % 
        % Input:
        %    mode    - up: iterates from leaves to root        
        %              down: iterates from root to leaves
        %              interior_up: iterates from base of tree (ignoring
        %              leaves) to root
        %              interior_down: iterates from root to base of tree,
        %              ignoring leaves
        %              leaves: iterates over leaves only
        % 
        %   T        - cell array of structs conforming to dimTree. Used to
        %              store/retrieve values associated to each node (e.g.
        %              the result of dimTree.emptyTree())
           if ~ischar(mode)
              error('mode must be a valid string'); 
           end
           if nargin == 2
               T = [];
           end
           switch mode
             case 'up'
               itr = dimTreeItrUp(dimTree,T);
             case 'down'
               itr = dimTreeItrDown(dimTree,T);
             case 'interior_up'
               itr = dimTreeItrIntUp(dimTree,T);
             case 'interior_down'
               itr = dimTreeItrIntDown(dimTree,T);
             case 'leaves'
               itr = dimTreeItrLeaves(dimTree,T);
           end
        end
        
        function T = copyCell(dimTree, T, S, copyName)
        % Copies the values of a cell array in to another cell array with 
        % the same tree structure
        %
        % Use:
        %    T = dimTree.copyCell(T,S,copyname)
        %
        % Input:
        %    T - target tree with the same structure as dimTree
        %    S - source tree (with the same structure as T)
        %    copyName - name of the variable to assign in T
        %
        % Output:
        %    T{i}{j} satisfies
        %    
        %     T{i}{j}.{copyName} = S{i}{j}
        %
        %  for each (i,j) index in dimTree
            itr = dimTree.iterator('down',T);
            while itr.advance()
                itr.setValue(copyName,S{itr.depth}{itr.idx}); 
            end
            T = itr.T;
        end
        
        function T = copyTree(dimTree, T, S, label_func)
        % Conditionally copies the values of a source cell array in to the target
        % cell array with the same tree structure.         
        % 
        % Use:
        %    T = dimTree.copyTree(T,S,{label_func})
        %
        % Input:
        %    T - target cell array
        %    S - source cell array
        %    label_func - single-input boolean function - returns true if
        %    the provided label should be copied (default: @(x) true)
        %
        % Output:
        %    T - target cell array with copied variables
        
           if exist('label_func','var')==0
               label_func = @(x) true;
           end
           itr = dimTree.iterator('down',S);
           while itr.advance()
               labels = itr.getLabels();
               
               for i=1:length(labels)
                   if label_func(labels{i})
                       T{itr.depth}{itr.idx}.(labels{i}) = itr.getValue(labels{i});
                   end
               end               
           end
        end
        
        function [U,B] = extractParams(dimTree,T, Uname, Bname)
        % Extracts parameters associated to each node of a cell array.
        %
        % [U,B] = dimTree.extractParams(T,{Uname},{Bname})
        %
        % Input: 
        %    T     - cell array of structs, same structure as dimTree
        %    Uname - variable name string associated to the leaves (default:
        %    'U')
        %    Bname - variable name string associated to the interior
        %    nodes (default: 'B')
        %
        % Output:
        %   [U,B] - collected parameters
            if exist('Uname','var')==0
                Uname = 'U';
            end
            if exist('Bname','var')==0
                Bname = 'B';
            end
            U = cell(dimTree.d,1);
            B = cell_skeleton(T);
            B = B{1:end-1};
            for i=1:length(U)
                I = dimTree.leaf_indices(:,i);
                U{i} = T{I(1)}{I(2)}.(Uname);
            end
            for i=1:dimTree.depth-1
                for j=1:length(T{i})
                    if ~dimTree.isLeaf(i,j)
                        B{i}{j} = T{i}{j}.(Bname);
                    end
                end
            end
        end
        
        function G = extractParam(dimTree,T,name)
        % Extracts a cell array tree structure from T containing one
        % variable from each node.
        %
        % Usage:
        %   G = dimTree.extractParam(T,name);
        %
        % Input:
        %   T - cell array of structs corresponding to dimTree
        %   name - variable name to extract
        %
        % Output:
        %   G - cell array whos entries G{i}{j} are the extracted variable
        %       values
        
           G = cell_skeleton(T);
           for i=1:length(T)
               for j=1:length(T{i})
                   G{i}{j} = T{i}{j}.(name);
               end
           end
        end
        
        function [U,B] = truncate(dimTree,X,rSVD)
        % Truncates a full tensor to HT format using the ranks specified by
        % the dimension tree. The error bounds of the truncated tensor,
        % X_trunc, satisfy
        %
        %    |X - X_{trunc}|_2 \le sqrt(2d - 3) * |X - X_{best}|_2
        %
        % Use:
        %  x = dimTree.truncate(X,{rSVD});
        %  [U,B] = dimTree.truncate(X,{rSVD});
        %
        % Input:
        %   X    - input tensor (must have compatible dimensions with dimTree)
        %   rSVD - true if randomized SVDs should be used, for
        %          large tensors X (default: false)
        % Output:
        %   x - vectorized parameters corresponding to the truncated tensor
        %   [U,B] = unvectorized parameters

            if exist('rSVD','var')==0
               rSVD = false; 
            end
            [U,B] = dimTree.fromVec(dimTree.randn());
            T = dimTree.emptyTree(U,B);
            distMode = isdistributed(X) | iscodistributed(X);
            
            itr = dimTree.iterator('leaves',T);
            C_lower = X;
            while itr.advance()
                U = itr.getValue('U');
                [nt,kt] = size(U);
                dim = itr.getDims();
                [leftSV,~,~] = svd(matricize(X,dim) * matricize(X,dim)');
                U = leftSV(:,1:kt); itr.setValue('U',U);                
                C_lower = ttm(C_lower, U',dim);
                clear leftSV;
            end                        
            C_level = C_lower;                        
            dimPartition = {};
            d = 1;
            for i=1:length(B{end})
                height = length(B);               
                if dimTree.isLeaf(height,i)
                    dimPartition{end+1} = [d];
                    d = d + 1;
                else
                    dimPartition{end+1} = [d,d+1];
                    d = d + 2;
                end
            end
            dimOffset = 1;
            itr = dimTree.iterator('up',itr.T);
            while itr.advance()
                if ~itr.isRoot()
                    if ~itr.isLeaf()                        
                        B = itr.getValue('B');
                        kt = size(B,3);
                        if rSVD
                            [leftSV,~,~] = randsvd(matricize(C_level,dimPartition{dimOffset}),kt+10);
                        else
                            [leftSV,~,~] = svd(matricize(C_level,dimPartition{dimOffset}));
                        end
                        B = dematricize(leftSV(:,1:kt),size(B),[1 2]);
                        itr.setValue('B',B);
                        clear leftSV;
                        beforeDims = 1:dimOffset-1;
                        afterDims = dimOffset+2:length(size(C_lower));                        
                        
                        C_lower = ttt(B,C_lower,[1 2],dimOffset:dimOffset+1);
                        % Consistent dimension ordering
                        if ~isempty(afterDims)
                            afterDims = afterDims-1;
                        end
                        if ~isempty(beforeDims)
                            beforeDims = beforeDims+1; 
                        end
                        C_lower = permute(C_lower,[beforeDims 1 afterDims]);
                        
                    end
                    dimOffset = dimOffset + 1;
                    if itr.idx == length(dimTree.nodes{itr.depth})
                        level = itr.depth;
                        C_level = C_lower;
                        dimOffset = 1;
                        if level > 2
                            dimPartition = {};
                            d = 1;
                            for i=1:length(dimTree.nodes{level-1})
                                if dimTree.isLeaf(level-1,i)
                                    dimPartition{end+1} = [d];
                                    d = d + 1;
                                else
                                    dimPartition{end+1} = [d,d+1];
                                    d = d + 2;
                                end
                            end
                        end
                    end
                else
                   B = C_level; 
                   itr.setValue('B',B);
                end                
            end
            [U,B] = dimTree.extractParams(itr.T);
            if nargout == 1
               U = dimTree.toVec(U,B); 
            end
        end        
        
        function G = gramian(dimTree,x,isOrthog)
        % Computes the Gramian matrices for an hierarchical tucker
        % tensor given from parameter matrices x
        % The Gramian matrix associated to a non-root node t
        % satisfies
        %    \lambda(G_t)_i = \sigma(X^(t))_i^2
        % where \lambda(G_t)_i is the ith eigenvalue of G_t 
        % and \sigma(X^(t))_i is the ith singular value of the
        % matricization X^(t)
        %
        % Use:
        %   G = dimTree.gramian(x)
        % Input:
        %   x       - vectorized parameter matrices (orthogonalized)
        % 
        % Output:
        %   G       - cell array with the same structure as
        %             dimTree, with each entry containing the Gramian associated to it            
                
            [U,B] = dimTree.fromVec(x);           
            if exist('isOrthog','var')==0
                isOrthog = true;
            end
            if ~isOrthog               
                error('Not implemented');
            end
            T = dimTree.emptyTree(U,B);
            itr = dimTree.iterator('interior_down',T);            
            while itr.advance()
               if itr.isRoot()
                   G = 1;
                   itr.setValue('G',1);
               else
                   G = itr.getValue('G');
               end
               B = itr.getValue('B');
               B_mod = ttm(B,G,3);
               itr.setLeftValue('G',ttt(B,B_mod,[2 3],[2 3]));
               itr.setRightValue('G',ttt(B,B_mod,[1 3],[1 3]));               
            end            
            G = dimTree.extractParam(itr.T,'G');
            
        end       
        
        function x = toVecGram(dimTree,G)
        % Vectorizes a cell array of Gramian matrices
            x = treeCell2Vec(G);
        end
               
        function G = fromVecGram(dimTree,x)
        % Converts a vectorized cell array of Gramian matrices back
        % to its original cell array form 
            G = cell_skeleton(dimTree.nodes);
            G{1}{1} = [1 1];
            for i=2:length(G)
                nodes = dimTree.nodes{i};
                for j=1:length(nodes)
                    k = nodes{j}.rank;
                    G{i}{j} = [k k];
                end
            end
            G = vec2TreeCell(x,G);
        end                
       
        function y = rank(dimTree,i,j,rank)
        % Function to either set or get the rank associated to a
        % particular node 
        %
        % Usage:
        %   a) rank = dimTree.rank(i,j)
        %   b) dimTree.rank(i,j,rank)
        %
        % Input: 
        %   i         - depth of the target node in the tree
        %   j         - horizontal index of the target node at
        %               depth level i
        %(b) rank     - new rank of the node at (i,j)
        %  
        %  Output:
        %(a) rank     - current rank of the node at (i,j)
            if nargin == 3
                %Get rank
                y = dimTree.nodes{i}{j}.rank;
            elseif nargin == 4
                if i==1 && j==1
                    error('root node rank must be 1');
                end
                dimTree.nodes{i}{j}.rank = rank;                                  
            else
               error('Invalid number of arguments'); 
            end
        end
        
        function z = isEmpty(dimTree,i,j)        
           z = isempty(dimTree.nodes{i}{j}); 
       end
       
        function z = isLeaf(dimTree,i,j)
        % Returns true if the index pair (i,j) corresponds to a leaf node,
        % false otherwise.
            z = length(dimTree.nodes{i}{j}.dims)==1;
        end
        function v = toVec(dimTree,U,B)
        % Vectorizes a set of parameter matrices in cell array form. 
        % 
        %   v = dimTree.toVec(U,B);
        % 
        % Input:
        %   U - Cell array of leaf matrices
        %   B - Cell array of transfer tensors
        % 
        % Output:
        %   v - vectorized parameters
            len = dimTree.sizeOfVec();
            v = zeros(len,1);
            total = 1;
            for i=1:length(U)
                len = numel(U{i});
                v(total:total+len-1) = vec(U{i});
                total = total + len;
            end
            for i=1:length(B)
                for j=1:length(B{i})
                    if ~isempty(B{i}{j})
                        len = numel(B{i}{j});
                        v(total:total+len-1) = vec(B{i}{j});
                        total = total + len;
                    end
                end
            end            
        end                                              

        function len = sizeOfVec(dimTree)
        % Size of the vectorized cell array of leaf matrices +
        % transfer tensors            
            len = dimTree.sizeOfUVec() + dimTree.sizeOfBVec();
        end
        
        function v = randn(dimTree)
        % Returns an (unorthogonalized) set of random hierarchical tucker parameter
        % matrices (in vector form)
        % 
        % Use:
        %   x = dimTree.randn()
        %
        % Output:
        %   x - random vectorized HT parameters
            len = dimTree.sizeOfVec();
            v = randn(len,1);            
            [U,B] = dimTree.fromVec(v);
            B{1}{1} = B{1}{1}/norm(B{1}{1});
            v = dimTree.toVec(U,B);
        end
        
        function v = randnc(dimTree)
            len = dimTree.sizeOfVec();
            v = randn(len,1) + 1i*randn(len,1);
        end
        
        function X = fullTree(dimTree,U,B, computeRoot, distributed)
        % Computes the full matrices associated with each node
        % 
        % Usage:
        %   X = dimTree.fullTree(x);
        %   X = dimTree.fullTree(U,B,{computeRoot}, {distributed});
        %
        % Input:
        %   x            - vectorized parameter matrices
        %   U            - leaf matrices (in cell array form)
        %   B            - transfer tensors (in cell array form)
        %   computeRoot  - true: compute the full tensor
        %                  (associated to the root node) (default: false)
        %   distributed  - true: if computeRoot is true, the root tensor is
        %                  distributed in parallel
        %
        %  Output:
        %   X            - cell array containing the (full) U matrices and
        %                  B transfer tensors associated to each node        
            if(exist('computeRoot','var')==0)
                computeRoot = false;
            end           
            if(nargin == 2)
               [U,B] = dimTree.fromVec(U); 
            end                       
            T = dimTree.emptyTree(U,B);            
            itr = dimTree.iterator('interior_up',T);
            while itr.advance()
               Ul = itr.getLeftValue('U');
               Ur = itr.getRightValue('U');
               B = itr.getValue('B');
               if itr.isRoot()
                   if computeRoot
                       if exist('distributed','var') && distributed
                           U = dimensionTree.construct_node(Ul,Ur,B,true,false);
                       else
                           U = dimensionTree.construct_node(Ul,Ur,B,false,false);
                       end
                   else
                       U = [];
                   end
               else
                   U = dimensionTree.construct_node(Ul,Ur,B,false,false);
               end
               clear Ul Ur B;
               itr.setValue( 'U' , U );                
            end
            X = itr.T;
            
        end                
        
        function X = full(dimTree, U,B,distributed,implicit)
        % Computes the full tensor from parameter matrices
        % 
        % Usage:
        %   X = dimTree.full(x);
        %   X = dimTree.full(U,B,{distributed});
        %
        % Input:
        %   x            - vectorized parameter matrices
        %   U            - leaf matrices (in cell array form)
        %   B            - transfer tensors (in cell array form)
        %   distributed  - true: full tensor is distributed (default: false)
        %
        %  Output:
        %   X            - full tensor
            if nargin == 1
               error('Not enough inputs'); 
            end
            if(nargin == 2)
               [U,B] = dimTree.fromVec(U); 
            end                       
            if exist('distributed','var')==0
                distributed = false;
            end
            if exist('implicit','var')==0
                implicit = false;
            end
            T = dimTree.emptyTree(U,B);            
            itr = dimTree.iterator('interior_up',T);
            while itr.advance()
               Ul = itr.getLeftValue('U');
               Ur = itr.getRightValue('U');
               B = itr.getValue('B');
               if itr.isRoot() 
                   U = dimensionTree.construct_node(Ul,Ur,B,distributed,false);
               else 
                   U = dimensionTree.construct_node(Ul,Ur,B,false,implicit);
               end
               clear Ul Ur B;
               itr.setLeftValue('U',[]);
               itr.setRightValue('U',[]);
               itr.setValue( 'U', U );    
            end
            X = itr.T{1}{1}.U;            
        end
        
        function X = fullND(dimTree,x,distributed,implicit)
        % Convenience function for expanding a full HT tensor + reshaping it
        %
        %   X = dimTree.fullND(x);
        %   X = dimTree.fullND(U,B,{distributed});
        %
        % are the same as, respectively,
        %   
        %   X = reshape(dimTree.full(x),dimTree.dims);
        %   X = reshape(dimTree.full(U,B,{distributed}),dimTree.dims);        
       
           if nargin == 1
              error('Not enough inputs'); 
           end
           [U,B] = dimTree.fromVec(x); 
           if exist('distributed','var')==0
               distributed = false; 
           end
           if exist('implicit','var')==0
               implicit = false;
           end
           X = dimTree.full(U,B,distributed,implicit);
           if distributed
               X = pSPOT.utils.distVec2
           else
               X = reshape(X,dimTree.dims);
           end
       end
       
        function X = fullDist(dimTree,x,implicit)
        % Convenience function for expanding a full HT tensor in a distributed manner        
        %   X = dimTree.fullDist(x);
        % is the same as
        %   [U,B] = dimTree.fromVec(x);
        %   X = dimTree.full(U,B,true);
            if exist('implicit','var')==0
                implicit = false;
            end
            [U,B] = dimTree.fromVec(x);
            X = dimTree.full(U,B,true,implicit);
        end        
        
    end
    
    methods (Access = protected)                
         function len = sizeOfUVec(dimTree)
         % Returns the length of the vector of U parameter matrices
             len = 0;
             for i=1:dimTree.d
                 idx = dimTree.leaf_indices(:,i);
                 
                 len = len + (dimTree.nodes{idx(1)}{idx(2)}.rank * dimTree.dims(i));
             end
         end
         function len = sizeOfBVec(dimTree)
         % Size of the vectorized cell array of transfer tensors
             len = 0;
             for i=dimTree.depth-1:-1:1                
                 for j=1:length(dimTree.nodes{i})                    
                     if ~dimTree.isLeaf(i,j)
                         left_idx = 2*(j-1)+1; right_idx = 2*(j-1)+2;
                         dims = [dimTree.nodes{i+1}{left_idx}.rank,...
                                 dimTree.nodes{i+1}{right_idx}.rank,...
                                 dimTree.nodes{i}{j}.rank];
                        
                         len = len + prod(dims);
                     end
                 end
             end
         end
    end    
    
    
    methods(Static)
        function Z = construct_node(U_left, U_right, B, distribute, implicit)
            if(exist('distribute','var')==0)
                distribute = false;
            end                        
            [kleft, kright, k] = size(B);  
            nl = size(U_left,1); nr = size(U_right,1);
            if(~distribute)
                if implicit && k > 1                                        
                    Z = opKron(U_right,U_left) * opMatrix(matricize(B,[1 2]));
                else
                    Z = opKron(U_right,U_left) * matricize(B,[1 2]);
                end
                %Z = matricize(ttm(B, {U_left,U_right},[1 2]),[1,2]);
            else
                P = oppKron2Lo(U_right, U_left);
                Z = P * vec(B);                
            end
         end
       
    end
    
end