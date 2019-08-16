classdef oppStack < oppSpot
    %OPPSTACK  Stack of vertically concatenated operators in parallel
    %
    %   S = oppStack([WEIGHTS], OP1, OP2, ...OPn,GATHER) creates a stacked 
    %   operator A consisting of the vertical concatenation of all 
    %   operators. When applied the operators are divided amongst the labs 
    %   and applied locally on each lab.
    %
    %   S = oppStack(N,OP) creates a stacked operator A using N number of
    %   repeating operators OP.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results.
    %
    %   Optional WEIGHTS vector:
    %
    %               [WEIGHT1*OP1
    %                WEIGHT2*OP2
    %                   ...
    %                WEIGHTn*OPn]
    %
    %   If the same weight is to be applied to each operator, set
    %   WEIGHTS to a scalar. When WEIGHTS is empty [], it is set to
    %   one. The WEIGHT parameter can be omitted as long as OP1 is a
    %   Spot operator; if not there is no way to
    %   decide whether it is a weight vector or operator.
    %
    %   *Note - only spot operators can be used in an oppStack, pSpot
    %   operators act in parallel already, and cannot be used with
    %   oppDictionary.
    %
    %   **Note - As of now the operators will be distributed according to
    %   the Matlab default codistribution scheme.
    %   for more info, type 'help codistributor1d'
    %
    %   See also oppBlockDiag, oppNumBlockDiag, oppDictionary
    
    %   Nameet Kumar - Oct 2010
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppStack(varargin)
            
            % Check parallel pool
            assert(parpool_size() > 0, 'parallel pool is not on');
            
            % Check for gather parameter
            if isscalar( varargin{end} ) && ~isa(varargin{end},'opSpot')
                gather = varargin{end};
                varargin(end) = [];
            else
                gather = 0;
            end
            
            % Check for weights
            nargs = length(varargin);
            
            if isnumeric(varargin{1}) % weights
                weights = varargin{1};
                weights = weights(:);
                
                if nargs == 2 % Repeating ops                    
                    if spot.utils.isposintscalar(varargin{1}) % repeat N times
                        weights = ones(weights,1);                        
                    end % Else: Repeating as many times as there are weights
                    
                    for i = 3:length(weights)+1
                        varargin{i} = varargin{2};
                    end
                    
                else % Non-repeating ops                    
                    if isscalar(varargin{1}) % Same weight applied to all
                        weights = weights*ones(nargs-1,1);
                        
                    else
                        assert(length(varargin{1}) == nargs-1,...
                            'Weights size mismatch');% Wrong weight size
                        % Else: Normal weights with normal ops
                    end
                end
                varargin(1) = []; % delete weights
                
            else    % no weights
                % Check for empty children
                nargs   = sum(~cellfun(@isempty,varargin));
                weights = ones(nargs,1);
            end
            
            % Standard pSpot checking and setup sizes
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            assert( all(n == n(1)), 'Operator sizes are not consistant');
            
            % Construct
            op = op@oppSpot('pStack', sum(m), n(1));
            op.cflag       = cflag;
            op.linear      = linear;
            op.children    = pSPOT.utils.compositeDef(opList);
            op.weights     = pSPOT.utils.compositeDef(weights);
            op.sweepflag   = true;
            op.gather      = gather;
            op.precedence  = 1;
            op.opsm        = m;
            
        end %Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Initialize
            opchildren = [op.children{:}];
            str = ['[',char(opchildren{1})];
            
            for ops=opchildren(2:end)
                str = strcat(str, '; ', char(ops{1}));
            end
            
            str = [str, ']'];
            
        end % Display
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = double(op)
            %OPPSTACK.DOUBLE Distributed doubling of oppStack
            %    A = double(op) will apply double to each child operator
            %    in oppStack, and return a distributed stack of explicit
            %    operators.
            
            loc_childs = op.children;
            loc_wgts   = op.weights;
            y_size     = [op.m op.n];
            spmd
                % Setup distribution stuffs
                chi_m = 0;
                y_part = codistributed.zeros(1,numlabs);
                
                % Preallocate
                y      = zeros(chi_m,y_size(2));
                if ~isempty(loc_childs)
                    y = double(opStack(loc_wgts,loc_childs{:}));
                    y_part(labindex) = chi_m;
                end
                y_part = gather(y_part);
                y_cod  = codistributor1d(1,y_part,y_size);
                y      = codistributed.build(y,y_cod,'noCommunication');
            end % spmd            
        end % double
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = drandn(A,Ncols)
            ncols = 1;
            if nargin == 2 % for easy multivectoring
                ncols = Ncols;
            end
            
            n = A.n;
            if isreal(A)
                x = randn(n,ncols);
            else
                x = randn(n,ncols) + 1i*randn(n,ncols);
            end
            
        end % drandn
        
    end % Methods    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if mode == 1
                % Checking distribution of x
                if isdistributed(x)
                    spmd, x_cod = getCodistributor(x); end
                    x_cod = x_cod{1};
                    assert(x_cod.Dimension ~= 1,...
                        'x cannot be distributed along first dimension');
                    x = gather(x);
                end

                % Mode 1
                % Setting up class variables and partition sizes
                loc_childs = op.children;
                loc_wgts   = op.weights;
                y_size     = [op.m size(x,2)]; % final global size

                spmd
                    % setup y parts
                    y_part = codistributed.zeros(1,numlabs);
                    
                    if ~isempty(loc_childs)                    
                        % Multiply
                        for i=1:length(loc_childs)
                            loc_childs{i} = loc_wgts(i) * loc_childs{i};
                        end
                        y = opStack(loc_childs{:})*x;
                    else
                        y = zeros(0,y_size(2), class(x));
                    end
                    
                    % Fill in the yparts
                    y_part(labindex) = size(y,1);

                    % Check for sparsity
                    aresparse = codistributed.zeros(1,numlabs);
                    aresparse(labindex) = issparse(y);                
                    if any(aresparse), y = sparse(y); end;

                    % Concatenating the results and distribute
                    y_part = gather(y_part);
                    y_cod  = codistributor1d(1,y_part,y_size);
                    y = codistributed.build(y,y_cod,'noCommunication');
                end % spmd mode 1
                
            else % mode == 2 (copied from oppDictionary)
                % X must be distributed
                assert(isdistributed(x),'X must be distributed');

                % Checking size of x
                spmd, x_cod = getCodistributor(x); end
                x_cod    = x_cod{1};
                x_part   = x_cod.Partition;
                chi_part = pSPOT.utils.defaultDistribution(length(op.opsm));
                assert(x_cod.Dimension == 1,... % Dimensional check
                    'x is not distributed along dimension 1');

                % Size checkings
                chi_num = 0;
                for i=1:parpool_size() 
                    chi_m = sum(op.opsm(chi_num+1:(chi_num+chi_part(i))));
                    assert(chi_m == x_part(i),...
                        'x size mismatch at lab %d, check your distribution',i);
                    chi_num = chi_num + chi_part(i);
                end
                loc_childs = op.children;
                loc_wgts   = op.weights;

                % Mode 1
                % Setting up preallocation size
                y_size = [op.n size(x,2)];

                spmd
                    % Setting up local parts
                    loc_x = getLocalPart(x);

                    % Preallocate y
                    y = zeros(y_size, class(loc_x));

                    if ~isempty(loc_childs)
                        for i=1:length(loc_childs)
                            % conjugate transpose
                            loc_childs{i} = conj(loc_wgts(i)) *...
                                            loc_childs{i}';
                        end
                        y = opDictionary(loc_childs{:}) * loc_x;
                    end

                    % Check for sparsity
                    aresparse            = codistributed.zeros(1,numlabs);
                    aresparse(labindex)  = issparse(y);
                    if any(aresparse), y = sparse(y); end;

                    % Summing the results and distribute
                    y = gplus(y,1); % The result is now on lab 1
                    y = codistributed(y,1,codistributor1d(2));

                end % spmd mode 2
            end
            
            if op.gather
                y = gather(y);
            end    %if we gathered, the data is on master client            
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Sweepable
            x = matldivide(op,b,mode);
        end % divide
    end % Protected Methods    
end % Classdef
