classdef oppSweep < oppSpot
    %OPPSWEEP   Operator that sweeps across parallel multivectors/arrays
    %
    %   B = oppSweep(OP,GATHER) allows operator OP to perform 
    %   sweeping multiplication across vectors or arrays distributed across
    %   the last dimension (which it is distributed), in parallel. 
    %   OP cannot be a pSpot operator.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppSweep(varargin)
                        
            % Check for parallel pool
            if parpool_size() == 0
                error('Parallel pool is not on');
            end
            
            % Settin' up the variables
            gather = 0;
            
            % Extract gather
            if isscalar(varargin{end}) && any(varargin{end} == [0 1 2 3])
                gather = varargin{end};
                varargin(end) = [];
            end
                        
            % Check for number of operators
            assert(length(varargin) == 1, 'Only one operator is supported');
            
            % Standard checking
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{1});
            
            % Construct
            op = op@oppSpot('pSweep', m, n);
            op.cflag     = cflag;
            op.linear    = linear;
            op.children  = opList;
            op.sweepflag = true;
            op.gather    = gather;
            
        end % constructor
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            
            str = ['pSweep(',char(op.children{1}),')'];
            
        end % Display
        
    end % methods
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Check for sizes
            if mode == 1
                if op.n ~= size(x,1)
                    error(...
               'Matrix dimensions must agree when multiplying by %s.',...
               char(op));
                end
            else
                if op.m ~= size(x,1)
                    error(...
               'Matrix dimensions must agree when multiplying by %s.',...
               char(op));
                end
            end
            
            % Setup variables
            A = op.children{1};
            opm = op.m; opn = op.n;
            nlabs = parpool_size();
            size_x = size(x);
            size_rest = size_x(2:end);
            
            spmd
                % Setup local parts
                local_x = getLocalPart(x);
                
                % Setup final codistributor
                finpart = codistributed.zeros(1,nlabs);
                if mode == 1
                    fingsize = [opm size_rest];
                else
                    fingsize = [opn size_rest];
                end
                
                % Multiply
                % Setup partition
                finpart(labindex) = size(local_x,length(size_x));
                
                % Vec local_x
                vec_x = local_x(:);
                % Setup kron
                size_loc = size(local_x);
                if mode == 1
                    kronop = opKron(opDirac(prod(size_loc(2:end))),A);
                else
                    kronop = opKron(opDirac(prod(size_loc(2:end))),A');
                end
                
                % Multiply
                local_y = kronop*vec_x;
                
                % Reshape
                size_y = size(local_x);
                if mode == 1
                    size_y(1) = opm;
                else
                    size_y(1) = opn;
                end
                
                local_y = reshape(local_y,size_y);
                
                % Check for sparsity
                aresparse = codistributed.zeros(1,numlabs);
                aresparse(labindex) = issparse(local_y);
                % labBarrier;
                if any(aresparse), local_y = sparse(local_y); end;
                
                fincodist = codistributor1d(length(size_x),finpart,fingsize);
                
                % Build codistributed y
                y = codistributed.build(local_y,fincodist,'noCommunication');
            end
            
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
        function x = divide(op,b,mode)
            % Sweepable
            x = matldivide(op,b,mode);
        end % divide
        
    end % Protected methods
    
end % classdef


















