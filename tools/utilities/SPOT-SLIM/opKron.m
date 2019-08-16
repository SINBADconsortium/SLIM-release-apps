classdef opKron < opSpot
%opKron   Kronecker tensor product.
%
%   opKron(OP1,OP2,...OPn) creates an operator that is the Kronecker
%   tensor product of OP1, OP2, ..., OPn.

%   Copyright 2009, Rayan Saab, Ewout van den Berg and 
%   Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        permutation; %Permutation vector of intergers defining the order to
        %use when the operators (children) of the Kronecker product are
        %applied to a data vector.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opKron(varargin)
            narg=nargin;
            
            %Test the case where varargin is a list
            if narg == 1
                narg=length(varargin{1});
                varargin=varargin{1};
            end
            
            if narg < 2
                error('At least two operators must be specified')
            end
            
            % Input matrices are immediately cast to opMatrix.
            for i=1:narg
                if isa(varargin{i},'numeric')
                    varargin{i} = opMatrix(varargin{i});
                elseif ~isa(varargin{i},'opSpot')
                    error('One of the operators is not a valid input.')
                end
            end
            
            % Determine operator size and complexity (this code is
            % general for any number of operators)
            opA       = varargin{1};
            [m,n]     = size(opA);
            cflag     = opA.cflag;
            linear    = opA.linear;
            sweepflag = opA.sweepflag;
            
            for i=2:narg
                opA       = varargin{i};
                cflag     = cflag  | opA.cflag;
                linear    = linear & opA.linear;
                sweepflag = sweepflag & opA.sweepflag;
                [mi,ni]   = size(opA);
                m = m * mi; n = n * ni;
            end
            
            % Construct operator
            op = op@opSpot('Kron', m, n);
            op.cflag       = cflag;
            op.linear      = linear;
            op.sweepflag   = sweepflag;
            op.children    = varargin;
            op.permutation =(1:narg);
            
            %Evaluate the best permutation to use when a multiplication is
            %applied
            if ~ (m == 0 || n == 0)
                op.permutation=op.best_permutation();
            end
            
            % Setting up implicit dimensions of output vector
            % Flipped
            varargin = fliplr(varargin);
            len      = length(varargin);
            op.ms    = cell(1,len);
            op.ns    = cell(1,len);
            for u = 1:len
                child_op = varargin{u};
                if length(child_op.ms) > 1
                    op.ms{u} = [op.ms{u} child_op.ms(:)'];
                else
                    op.ms{u} = [op.ms{u} child_op.ms{:}];
                end
                if length(child_op.ns) > 1
                    op.ns{u} = [op.ns{u} child_op.ns(:)'];
                else
                    op.ns{u} = [op.ns{u} child_op.ns{:}];
                end
            end
                
        end % Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            str=['Kron(',char(op.children{1})];
            
            % Get operators
            for i=2:length(op.children)
                str=strcat(str,[', ',char(op.children{i})]);
            end
            str=strcat(str,')');
        end % Char
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = headerMod(op,header,mode)
            % Extract explicit size indices
            exsize = header.exsize;
            href   = @spot.data.headerRef; % Used function handles because
            hasg   = @spot.data.headerAsgn; % its shorter

            if mode == 1
                
                % Setup variables
                opList = fliplr(op.children); % Last op applied first
                % Number of output dimensions
                n_out_dims = length(spot.utils.uncell(op.ms)) +...
                                                        size(exsize,2) - 1;
                
                % Preallocate and setup header
                header_out        = header;
                header_out.dims   = n_out_dims;
                header_out.size   = zeros(1,n_out_dims);
                header_out.origin = zeros(1,n_out_dims);
                header_out.delta  = zeros(1,n_out_dims);
                header_out.unit   = zeros(1,n_out_dims);
                header_out.label  = zeros(1,n_out_dims);
                
                % Fill in the non-first-dimension header parts
                if size(exsize,2) > 1
                    h_part = href(header,exsize(1,2):exsize(end,end));
                    header_out = hasg(header_out,h_part,...
                        length(spot.utils.uncell(op.ms))+1:...
                        length(header_out.size));
                end
                
                % Replace old first (collapsed) dimensional sizes with
                % operator sizes.
                i = 1;
                x = 1;
                % Extract collapsed first dims from header.
                first_header = href(header,exsize(:,1));
                for u = 1:length(op.children)
                    % Input header (including collapsed)
                    y            = length(spot.utils.uncell(op.ns{u}))+x-1;
                    in_header    = href(first_header,x:y);
                    in_header.exsize = [1;y-x+1];
                    % child header
                    child_header = headerMod(opList{u},in_header,mode);
                    % Assignment indices
                    j            = length(spot.utils.uncell(op.ms{u}))+i-1;
                    % header assignment
                    oldsize      = length(header_out.size);
                    header_out   = hasg(header_out,child_header,i:j);
                    newsize      = length(header_out.size);
                    i            = j + 1 + newsize - oldsize;
                    x            = y + 1;
                end
                exsize_out = 1:length(header_out.size);
                exsize_out = [exsize_out;exsize_out];
                h          = header_out;
                h.exsize   = exsize_out;
            else
                % Setup variables
                opList = fliplr(op.children); % Last op applied first
                % Number of output dimensions
                n_out_dims = length(spot.utils.uncell(op.ns)) +...
                                                        size(exsize,2) - 1;
                
                % Preallocate and setup header
                header_out        = header;
                header_out.dims   = n_out_dims;
                header_out.size   = zeros(1,n_out_dims);
                header_out.origin = zeros(1,n_out_dims);
                header_out.delta  = zeros(1,n_out_dims);
                header_out.unit   = zeros(1,n_out_dims);
                header_out.label  = zeros(1,n_out_dims);
                
                % Fill in the non-first-dimension header parts
                if size(exsize,2) > 1
                    h_part = href(header,exsize(1,2):exsize(end,end));
                    header_out = hasg(header_out,h_part,...
                        length(spot.utils.uncell(op.ms))+1:...
                        length(header_out.size));
                end
                
                % Replace old first (collapsed) dimensional sizes with 
                % operator sizes.
                i = 1;
                x = 1;
                % Extract collapsed first dims from header.
                first_header = href(header,exsize(:,1));
                for u = 1:length(op.children)
                    % Input header (including collapsed)
                    y            = length(spot.utils.uncell(op.ms{u}))+x-1;
                    in_header    = href(first_header,x:y);
                    in_header.exsize = [1;y-x+1];
                    % child header
                    child_header = headerMod(opList{u},in_header,mode);
                    % Assignment indices
                    j            = length(spot.utils.uncell(op.ns{u}))+i-1;
                    % header assignment
                    oldsize      = length(header_out.size);
                    header_out   = hasg(header_out,child_header,i:j);
                    newsize      = length(header_out.size);
                    i            = j + 1 + newsize - oldsize;
                    x            = y + 1;
                end
                exsize_out = 1:length(header_out.size);
                exsize_out = [exsize_out;exsize_out];
                h          = header_out;
                h.exsize   = exsize_out;
            end
            
        end % headerMod
    end % Methods
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            % The Kronecker product (KP) is applied to the righthand matrix
            % taking in account the best order to apply the operators.
            % That necessitates to decompose the KP in successive matrix
            % products with terms of type I(a) kron A kron I(b).
            % A is the operator to apply. I(a) and I(b) are identity
            % matrices with respective sizes a and b.
            
            opList       = op.children; % Contains list of opKron children
            ncol         = size(x,2); % Number of columns of 'x'
            nbr_children = length(opList); % Number of children
            
            %%%%%%%%%%%%%%%%%%%%%%Multiplication%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mode == 1 % Classic mode
                ops = fliplr(opList);
                perm = op.permutation; % Permutation to take in account.
                perm = length(perm) - perm + 1;
                m = op.m; % Height of the resulting matrix
                
                dimOrder = 1:nbr_children;
                dimSize = zeros(1,length(ops));
                for ind = 1:length(ops)
                    dimSize(ind) = ops{ind}.n;
                end
                
                if ncol > 1
                    dimOrder = [dimOrder nbr_children+1];
                    dimSize(end+1) = ncol;
                end
                
                % first operator
                theOp = ops{perm(1)};
                % reshape into 2 groups
                if perm(1) == dimOrder(1) % if the first dim goes first
                    % x -> group 1: 1 ; group 2: 2:end
                    x = reshape(x,theOp.n,numel(x)/theOp.n);
                    dimOrder = {dimOrder(1) dimOrder(2:end)};
                    dimSize = {dimSize(1) dimSize(2:end)};
                else % perm(1) ~=dimOrder(1)
                    fetchIndex = find(dimOrder == perm(1));
                    % note: fetchIndex should have only one element
                    
                    % reshape:
                    %   group 1: 1:fetchIndex-1
                    %   group 2: fetchIndex:end
                    group1n = prod(dimSize(1:fetchIndex-1));
                    group2n = prod(dimSize(fetchIndex:length(dimOrder)));
                    
                    x = reshape(x,group1n,group2n);
                    
                    dimOrder = {dimOrder(1:fetchIndex-1) dimOrder(fetchIndex:length(dimOrder))};
                    dimSize = {dimSize(1:fetchIndex-1) dimSize(fetchIndex:length(dimSize))};
                    
                    x = x.';
                    x = reshape(x,theOp.n,numel(x)/theOp.n);
                    
                    dimOrder = {dimOrder{2}(1) [dimOrder{2}(2:end) dimOrder{1}]};
                    dimSize = {dimSize{2}(1) [dimSize{2}(2:end) dimSize{1}]};
                end
                
                for permInd = 2:length(perm)
                    if ~theOp.isDirac
                        x = applyMultiply(theOp, x, 1);
                        dimSize{1} = theOp.m;
                    end
                    theOp = ops{perm(permInd)};
                    if perm(permInd) == dimOrder{2}(1) % if the dimension is already the first element of the second group
                        x = x.';
                        x = reshape(x,theOp.n,numel(x)/theOp.n);
                        
                        dimOrder = {dimOrder{2}(1) [dimOrder{2}(2:end) dimOrder{1}]};
                        dimSize = {dimSize{2}(1) [dimSize{2}(2:end) dimSize{1}]};
                        
                    else % if not already the first element of the second group
                        % reshape:
                        %   group 1: 1:fetchIndex-1
                        %   group 2: fetchIndex:end
                        dimOrder = [dimOrder{1} dimOrder{2}];
                        dimSize = [dimSize{1} dimSize{2}];
                        fetchIndex = find(dimOrder == perm(permInd));
                        
                        group1n = prod(dimSize(1:fetchIndex-1));
                        group2n = prod(dimSize(fetchIndex:length(dimOrder)));
                        
                        x = reshape(x,group1n,group2n);
                        
                        dimOrder = {dimOrder(1:fetchIndex-1) dimOrder(fetchIndex:length(dimOrder))};
                        dimSize = {dimSize(1:fetchIndex-1) dimSize(fetchIndex:length(dimSize))};
                        
                        x = x.';
                        x = reshape(x,theOp.n,numel(x)/theOp.n);
                        
                        dimOrder = {dimOrder{2}(1) [dimOrder{2}(2:end) dimOrder{1}]};
                        dimSize = {dimSize{2}(1) [dimSize{2}(2:end) dimSize{1}]};
                    end
                end
                
                x = applyMultiply(theOp, x, 1);
                dimSize{1} = theOp.m;
                
                dimSize = [dimSize{1} dimSize{2}];
                
                dimOrder = [dimOrder{1} dimOrder{2}];
                
                lowerInd = find(dimOrder==1):length(dimOrder);
                higherInd = 1:find(dimOrder==1)-1;
                
                x = reshape(x,[prod(dimSize(higherInd)) prod(dimSize(lowerInd))]);
                x = x.';
                
                x = reshape(x,m,ncol);
            elseif mode == 2 % Transpose mode
                ops = fliplr(opList);
                perm = op.permutation; % Permutation to take in account.
                
                n = op.n; % Height of the resulting matrix
                
                dimOrder = 1:nbr_children;
                dimSize = zeros(1,length(ops));
                for ind = 1:length(ops)
                    dimSize(ind) = ops{ind}.m;
                end
                
                if ncol > 1
                    dimOrder = [dimOrder nbr_children+1];
                    dimSize(end+1) = ncol;
                end
                
                % first operator
                theOp = ops{perm(1)};
                % reshape into 2 groups
                if perm(1) == dimOrder(1) % if the first dim goes first
                    % x -> group 1: 1 ; group 2: 2:end
                    x = reshape(x,theOp.m,numel(x)/theOp.m);
                    dimOrder = {dimOrder(1) dimOrder(2:end)};
                    dimSize = {dimSize(1) dimSize(2:end)};
                else % perm(1) ~=dimOrder(1)
                    fetchIndex = find(dimOrder == perm(1));
                    % note: fetchIndex should have only one element
                    
                    % reshape:
                    %   group 1: 1:fetchIndex-1
                    %   group 2: fetchIndex:end
                    group1n = prod(dimSize(1:fetchIndex-1));
                    group2n = prod(dimSize(fetchIndex:length(dimOrder)));
                    
                    x = reshape(x,group1n,group2n);
                    dimOrder = {dimOrder(1:fetchIndex-1) dimOrder(fetchIndex:length(dimOrder))};
                    dimSize = {dimSize(1:fetchIndex-1) dimSize(fetchIndex:length(dimSize))};
                    
                    x = x.';
                    x = reshape(x,theOp.m,numel(x)/theOp.m);
                    dimOrder = {dimOrder{2}(1) [dimOrder{2}(2:end) dimOrder{1}]};
                    dimSize = {dimSize{2}(1) [dimSize{2}(2:end) dimSize{1}]};
                end
                
                for permInd = 2:length(perm)
                    if ~theOp.isDirac
                        x = applyMultiply(theOp, x, 2);
                        dimSize{1} = theOp.n;
                    end
                    
                    theOp = ops{perm(permInd)};
                    if perm(permInd) == dimOrder{2}(1) % if the dimension is already the first element of the second group
                        x = x.';
                        x = reshape(x,theOp.m,numel(x)/theOp.m);
                        dimOrder = {dimOrder{2}(1) [dimOrder{2}(2:end) dimOrder{1}]};
                        dimSize = {dimSize{2}(1) [dimSize{2}(2:end) dimSize{1}]};
                        
                    else % if not already the first element of the second group
                        % reshape:
                        %   group 1: 1:fetchIndex-1
                        %   group 2: fetchIndex:end
                        dimOrder = [dimOrder{1} dimOrder{2}];
                        dimSize = [dimSize{1} dimSize{2}];
                        fetchIndex = find(dimOrder == perm(permInd));
                        
                        group1n = prod(dimSize(1:fetchIndex-1));
                        group2n = prod(dimSize(fetchIndex:length(dimOrder)));
                        
                        x = reshape(x,group1n,group2n);
                        dimOrder = {dimOrder(1:fetchIndex-1) dimOrder(fetchIndex:length(dimOrder))};
                        dimSize = {dimSize(1:fetchIndex-1) dimSize(fetchIndex:length(dimSize))};
                        
                        x = x.';
                        x = reshape(x,theOp.m,numel(x)/theOp.m);
                        dimOrder = {dimOrder{2}(1) [dimOrder{2}(2:end) dimOrder{1}]};
                        dimSize = {dimSize{2}(1) [dimSize{2}(2:end) dimSize{1}]};
                    end
                end
                
                x = applyMultiply(theOp, x, 2);
                dimSize{1} = theOp.n;
                
                dimSize = [dimSize{1} dimSize{2}];
                
                dimOrder = [dimOrder{1} dimOrder{2}];
                
                lowerInd = find(dimOrder==1):length(dimOrder);
                higherInd = 1:find(dimOrder==1)-1;
                
                x = reshape(x,[prod(dimSize(higherInd)) prod(dimSize(lowerInd))]);
                x = x.';
                
                x = reshape(x,n,ncol);
            end
            
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,x,mode)
            % Depends on sweepflag
            if op.sweepflag
                x = matldivide(op,x,mode);
            else
                x = lsqrdivide(op,x,mode);
            end
        end % divide
        
    end %Protected methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % best_permutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Returns the best permutation associated to this Kronecker product
        function perm = best_permutation(op)
            list = op.children; % List of 'op''s children
            cost = zeros(length(list),2); % Computational costs of the
            % operators (children of 'op'). This is simply a numeric
            % representation of theirs shapes, which will affect 
            % computation time. Operators with low computational costs 
            % should be applied first.
            for i=1:length(list)
                % Cost = (nbr_rows-nbr_columns) / (size of the operator)
                cost(i,1) = (size(list{i},1)-size(list{i},2))/...
                    (size(list{i},1)*size(list{i},2));
                cost(i,2)=int8(i);
            end
            cost=sortrows(cost)';
            
            perm = cost(2,:);
        end
        
    end % private methods
end % opKron