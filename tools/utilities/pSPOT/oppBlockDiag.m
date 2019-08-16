classdef oppBlockDiag < oppSpot
    %OPPBLOCKDIAG   Operator-diagonal operator in parallel.
    %
    %   B = oppBlockDiag(OP1, OP2,...,OPN,GATHER) creates a compound
    %   block operator with the input operators OP1, OP2,... on the 
    %   diagonal of B, e.g., B = DIAG([OP1 OP2 ... OPN]). When multiplying 
    %   the operators are distributed among the labs and multiplied locally
    %   on each.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint 
    %            multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    %
    %   B = oppBlockDiag([WEIGHT],OP1,...,OPN,GATHER) additionally
    %   weights each block by the elements of the vector WEIGHT. If
    %   only a single operator is given it is replicated as many times
    %   as there are weights. If a scalar weight is given for multiple
    %   operators the scalar WEIGHT will be expanded to the size of the
    %   operators, ie. WEIGHT*ones(N,1)
    %
    %   B = oppBlockDiag(N,OP,GATHER) similar as above with WEIGHT
    %   equal to ones(N,1). This will cause operator OP to be repeated
    %   N times.
    %
    %   See also oppNumBlockDiag, oppDictionary, oppStack
    
    %   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
    %   See the file COPYING.txt for full copyright information.
    %   Use the command 'spot.gpl' to locate this file.
    
    %   http://www.cs.ubc.ca/labs/scl/spot
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppBlockDiag(varargin)
            
            % Check parallel pool
            assert(parpool_size() > 0, 'Parallel pool is not on');
            
            % Extract gather
            if isscalar(varargin{end}) && ~isa(varargin{end},'opSpot')
                gather        = varargin{end};
                varargin(end) = [];
            else
                gather        = 0;
            end
            
            % Extract weights and fill in repeating ops
            nargs = length(varargin);
            
            if isnumeric(varargin{1}) % weights
                weights = varargin{1};
                weights = weights(:);
                
                if nargs == 2 % Repeating ops
                    % repeating N times
                    if spot.utils.isposintscalar(varargin{1}) 
                        weights = ones(weights,1);                        
                    end % Else: Repeating as many times as there are weights
                    
                    for i = 3:length(weights)+1
                        varargin{i} = varargin{2};
                    end
                    
                else % Non-repeating ops                    
                    if isscalar(varargin{1}) % Same weight applied to all
                        weights = weights*ones(nargs-1,1);
                    else
                        if length(varargin{1}) ~= nargs-1
                            % Incorrect weight size
                            error('Weights size mismatch');
                        end
                        % Else: Normal weights with normal ops
                    end
                end
                varargin(1) = []; % delete weights
                
            else    % no weights
                % Check for empty children
                nargs   = sum(~cellfun(@isempty,varargin));
                weights = ones(nargs,1);
            end
            
            % Standard checking and setup sizes
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            
            % Construct operator
            op = op@oppSpot('pBlockDiag', sum(m), sum(n));
            op.cflag     = cflag;
            op.linear    = linear;
            op.children  = pSPOT.utils.compositeDef(opList);
            op.weights   = pSPOT.utils.compositeDef(weights);
            op.sweepflag = true;
            op.gather    = gather;
            op.opsn      = n;
            op.opsm      = m;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            loc_childs = [op.children{:}];
            str        = 'pBlockDiag(';
            
            if ~isnumeric( loc_childs{1} )
                for child = loc_childs
                    str = strcat(str,char(child{1}),', ');
                end
            else
                m = op.opsm; n = op.opsn;
                for i=1:length(loc_childs)
                    str = strcat(str,'Matrix(',int2str(m(i)),',', ...
                        int2str(n(i)),'), ');
                end
            end
            str = [str(1:end-2), ')'];
        end % Display
        
    end % Methods
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Checking x
            assert(isdistributed(x),'x is not distributed');
            
            % Checking size of x
            spmd, x_cod = getCodistributor(x); end
            x_cod    = x_cod{1};
            x_part   = x_cod.Partition;
            chi_part = pSPOT.utils.defaultDistribution(length(op.opsn));
            
            assert(x_cod.Dimension == 1,... % Dimensional check
                'x is not distributed along dimension 1');
            
            chi_num = 0;
            for i=1:parpool_size()
                chi_m = sum(op.opsm(chi_num+1:(chi_num+chi_part(i))));
                chi_n = sum(op.opsn(chi_num+1:(chi_num+chi_part(i))));
                
                if mode == 1
                    assert(chi_n == x_part(i),...
                        'x size mismatch at lab %d, check distribution',i);
                else % mode 2
                    assert(chi_m == x_part(i),...
                        'x size mismatch at lab %d, check distribution',i);
                end
                chi_num = chi_num + chi_part(i);
            end            
            
            % Setting up the variables and partition size
            loc_childs = op.children;
            loc_wgts   = op.weights;
            
            if mode ==  1
                y_size = [op.m size(x,2)]; % final global size
            else
                y_size = [op.n size(x,2)];
            end
            
            spmd
                % Setting up the local parts
                loc_x  = getLocalPart(x);
                
                % Setting up codistributor for y
                y_part = codistributed.zeros(1,numlabs);
                
                if ~isempty(loc_childs)                    
                    if length(loc_childs) == 1
                        % Just to avoid repeating op case
                        B = loc_childs{1};
                        if mode == 1
                            y = loc_wgts .* (B*loc_x);
                            y_part(labindex) = B.m;
                        else
                            y = conj(loc_wgts) .* (B'*loc_x);
                            y_part(labindex) = B.n;
                        end
                    else
                        B = opBlockDiag(loc_wgts,loc_childs{:});
                        if mode == 1
                            y                = B*loc_x;
                            y_part(labindex) = B.m;
                        else
                            y                = B'*loc_x;
                            y_part(labindex) = B.n;
                        end
                        
                    end
                else
                    y = zeros(0,size(x,2), class(loc_x));
                end
                
                % Check for sparsity
                aresparse            = codistributed.zeros(1,numlabs);
                aresparse(labindex)  = issparse(y);
                if any(aresparse), y = sparse(y); end;
                
                % Codistribute y
                y_cod = codistributor1d(1,y_part,y_size);
                y     = codistributed.build(y,y_cod,'noCommunication');
            end % spmd
            
            % Gather
            if mode == 1 
                if op.gather == 1 || op.gather == 2
                    y = gather(y);
                end
            else % mode == 2
                if op.gather == 1 || op.gather == 3
                    y = gather(y);
                end
            end % gather            
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = divide(op,x,mode)
            
            % Setting up the variables
            gather = op.gather;
            
            % Checking x
            assert(isdistributed(x),'x is not distributed');
            
            % Checking size of x
            spmd, x_cod = getCodistributor(x); end
            x_cod    = x_cod{1};
            x_part   = x_cod.Partition;
            chi_part = pSPOT.utils.defaultDistribution(length(op.opsn));
            n_labs   = parpool_size();
            
            assert(x_cod.Dimension == 1,... % Dimensional check
                'x is not distributed along dimension 1');
            
            chi_num = 0;
            for i=1:n_labs
                chi_m = sum(op.opsm(chi_num+1:(chi_num+chi_part(i))));
                chi_n = sum(op.opsn(chi_num+1:(chi_num+chi_part(i))));
                
                if mode == 1
                    assert(chi_m == x_part(i),...
                        'x size mismatch at lab %d, check distribution',i);
                else % mode 2
                    assert(chi_n == x_part(i),...
                        'x size mismatch at lab %d, check distribution',i);
                end
                chi_num = chi_num + chi_part(i);
            end
            
            % Setting up the variables and global sizes            
            loc_childs = op.children;            
            if mode ==  1
                y_size   = [op.n size(x,2)];
                loc_wgts = op.weights;
            else
                y_size   = [op.m size(x,2)];
                loc_wgts = pSPOT.utils.compositeDef(...
                                conj(vertcat(op.weights{:})));
            end
            
            spmd
                % Setting up the local parts
                loc_x = getLocalPart(x);
                
                % Setting up codistributor for y
                y_part = codistributed.zeros(1,numlabs);
                
                if ~isempty(loc_childs)                    
                    % Extracting local sizes
                    size_mn = cellfun(@size,loc_childs,'UniformOutput',false);
                    size_mn = [size_mn{:}];
                    loc_m   = sum(size_mn(1:2:end)); % sum odd
                    loc_n   = sum(size_mn(2:2:end)); % sum even
                    
                    if mode == 1
                        y_part(labindex) = loc_n;
                        
                        % Preallocate y
                        y = zeros(loc_n,size(x,2), class(loc_x));
                        
                        % Divide
                        j = 0; k = 0;
                        for i=1:length(loc_childs) % Divide by operator
                            child = loc_childs{i};
                            y(j+1:j+child.n,:) = loc_wgts(i) .* ...
                                          (child \ loc_x(k+1:k+child.m,:));
                            j     = j + child.n;
                            k     = k + child.m;
                        end
                        
                    else % mode 2
                        y_part(labindex) = loc_m;
                        
                        % Preallocate y
                        y = zeros(loc_m,size(x,2), class(loc_x));
                        
                        % Divide
                        j = 0; k = 0;
                        for i=1:length(loc_childs) % Divide by operator
                            child = ctranspose(loc_childs{i});
                            y(j+1:j+child.n,:) = loc_wgts(i) .* ...
                                (child \ loc_x(k+1:k+child.m,:));
                            j = j + child.n;
                            k = k + child.m;
                        end                        
                    end
                else
                    y = zeros(0,size(x,2));
                end
                
                % Check for sparsity
                aresparse           = codistributed.zeros(1,numlabs);
                aresparse(labindex) = issparse(y);
                if any(aresparse),y = sparse(y); end;
                
                % Codistribute y
                y_cod = codistributor1d(1,y_part,y_size);
                y     = codistributed.build(y,y_cod,'noCommunication');
                
            end % spmd
            
            % Gather
            if mode == 1
                if op.gather == 1 || op.gather == 2
                    y = gather(y);
                end
            else % mode == 2
                if op.gather == 1 || op.gather == 3
                    y = gather(y);
                end
            end % gather            
        end % divide
    end % Protected Methods    
end % Classdef
