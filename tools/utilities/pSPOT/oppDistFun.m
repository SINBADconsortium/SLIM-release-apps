classdef oppDistFun < oppSpot
    %OPPDISTFUN     Fun(ction) evaluation over slices of distributed arrays
    %
    %   Q = oppDistFun(A1,A2,...,AN,F,GATHER), where A1,A2,...,AN are 
    %   distributed arrays and F is a function handle that takes in local 
    %   parts of A1,A2,...,AN and gives a local part of the final answer. 
    %   The arguments of F has to be standardized as: 
    %   y = F(a1,a2,...,an,x,mode), where a1,a2,...,an corresponds to
    %   the local parts of A1,A2,...,AN, and x is the distributed vector 
    %   that the operator is applied on.
    %   mode = 1 defines the forward mode
    %   mode = 2 defines the adjoint mode
    %   mode = 0 will return the sizes, complexity and linearity in an
    %   array in the following format: [m n cflag linflag]
    %   m and n has to be the conceptual size of Q so that the
    %   multiplication sizes would match.
    %   cflag is the complexity of this operator.
    %   linflag is the linearity of this operator.
    %
    %   GATHER specifies whether to gather the results to a local array
    %   or leave them distributed, default is 0.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint 
    %            multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    %   
    %   ## Example use case:
    %   
    %   Computing expression over distributed index k
    %     y_k = (A_k * x_k - s_k)
    %   where
    %     A = size m by n by k
    %     x = size n by k
    %     s = size k
    %     (A, x, s are all distributed over the last dimension where size=k)
    % 
    %   We write a custom function F that performs (Ax-s) for each k
    %     F <-- function y = F(A, B, x, mode)
    %               if mode==1
    %               y = A * x - B;
    %               elseif mode==2
    %               y = A' * x - conj(B);
    %               elseif mode=='0'
    %               retrun {m, n}
    %           end
    % 
    %     P <---- constructed with `oppDistFun({A,s},F)`
    % 
    %     Then y = P*x will perform:
    % 
    %     for each A(:,k), s(k), k = 1:n (distributed)
    %        y(k) = F(A(:,k),s(k),x(k),1)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        fun;
        local_m;
        local_n;
        sizA;
        ndimA;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppDistFun(varargin)
                        
            if parpool_size() == 0 % Check for parallel pool
                error('Parallel pool is not open');
            end
            
            % Setup and extract variables
            opgather = 0;
            if isscalar(varargin{end}) && isnumeric(varargin{end})
                assert(any(varargin{end} == [0 1 2 3]),...
                    'Gather must be 0,1,2 or 3')
                opgather      = varargin{end};
                varargin(end) = [];
            end
            
            % Check for and extract function handle
            assert(isa(varargin{end},'function_handle'),...
                'F must be a function handle');
            F             = varargin{end};
            varargin(end) = [];
            
            % Turn off stupid warning
            warning('off','distcomp:codistributed:InvalidNumberOfLabs');
            
            % Store all ops as codistributed arrays inside cells            
            lenvar  = length(varargin);
            lastdim = zeros(1,lenvar);
            codist  = cell(1,lenvar);
            coddims = zeros(1,lenvar);
            codpart = cell(1,lenvar);
            ops = cell(1,lenvar);
            for i = 1:lenvar
               d = varargin{i};
               spmd 
                   ops{i} = d;
                   assert(iscodistributed(d), 'A must be distributed');
                   codd   = getCodistributor(d);
                   sized  = size(d);
               end
               sized      = sized{1};
               lastdim(i) = sized(end);
               codist{i}  = codd{1};
               coddims(i) = codist{i}.Dimension;
               if isempty(codist{i}.Partition)
                   codpart(i) = 0; % 0 means default partition
               else
                   codpart{i} = codist{i}.Partition;
               end
            end
            
            % Check for stuffs
            assert(all(lastdim == lastdim(1)),...
            'The last dimension must be of the same length');
            assert(all(coddims == ndims(varargin{1})),...
            'A must be distributed along the last dimension');
            for i=1:lenvar
                assert(all(codpart{i} == codpart{1}),...
            'Partition of distributed dimension must be the same');
            end
                                    
            % Extract parameters from function
            cell_padding = cell(1,nargin(F)-1);
            func_info    = F(cell_padding{:},0);
            m     = func_info(1); n       = func_info(2); 
            cflag = func_info(3); linflag = func_info(4);
            
            if ~spot.utils.isposintscalar(m) || ~spot.utils.isposintscalar(n) % check m and n
              error('Dimensions of operator must be positive integers.');
            end
            
            % Setup sizes
            sizA    = size(varargin{1});
            ndimA   = ndims(varargin{1});
            m_total = m*sizA(end);
            n_total = n*sizA(end);
            clear varargin;
                        
            % Construct oppDistFun
            op           = op@oppSpot('DistFun', m_total, n_total);
            op.local_m   = m;
            op.local_n   = n;
            op.children  = ops;
            op.fun       = F;
            op.cflag     = cflag;
            op.linear    = linflag;
            op.sweepflag = true;
            op.gather    = opgather;
            op.opsn      = n*ones(1,sizA(end));
            op.opsm      = m*ones(1,sizA(end));
            op.sizA      = sizA;
            op.ndimA     = ndimA;
            
        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
           str = 'oppDistFun';
            
        end % Display
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            if ~isa(op,'oppDistFun')
                error('Left multiplication not taken in account')
            else
                assert( isvector(x) , 'Please use vectorized matrix')
            end
            
            % Setup variables
            ops = op.children;
            F   = op.fun;
            if mode == 1
                xsize = op.local_n;
            else
                xsize = op.local_m;
            end
            ndimA = op.ndimA;
            
            % distribute x if necessary
            x = pSPOT.utils.scatterchk(x,mode,op.gather);
            % WARNING: do we check whether x and A1,A2,... are distributed equally?
            
            % Check for the distribution of x
            assert(isdistributed(x),'X must be distributed')
            
            spmd
                % Setup local parts
                xloc = getLocalPart(x);
                for i = 1:length(ops)
                   ops{i} = getLocalPart(ops{i});
                end
                
                % determine local size in the distributed dimension
                if ndims(ops{1}) < ndimA
                  % happens when loacal size == 1
                  nSliceA_local = 1;
                else
                  sizeA = size(ops{1});
                  nSliceA_local = sizeA(end);
                end
                
                if nSliceA_local > 0
                  % Setup y (not sure what the sizes of the other dims will be, so we initialize with cell arrays)
                  y = cell(1,nSliceA_local);
                
                  % Loop over the slices and apply F
                  n = 0;
                  for i=1:nSliceA_local
                      slice = cell(1,length(ops));
                      for j = 1:length(ops) % Get last-dimensional slice
                         slice{j} = pSPOT.utils.ndind(ops{j},ndimA,ndimA,i);
                      end
                    
                      % Get x slice
                      xslice = xloc(1 + n : xsize + n);
                      y{i} = F(slice{:},xslice,mode);
                      n = n + xsize;
                  end
                
                  % Stack y together
                  y = vertcat(y{:});
                else
                  % this worker recieved no slices, so will contribute an empty y (needs to be 2d because codistributor1d expects at least 2d)
                  y = zeros(0,1);
                end
                
                % Build y
                ypart = codistributed.zeros(1,numlabs);
                ypart(labindex) = size(y,1);
                ygsize = [sum(gather(ypart)) 1];
                ycod = codistributor1d(1,ypart,ygsize);
                y = codistributed.build(y,ycod,'noCommunication');
            end % spmd
            
            % gather
            y = pSPOT.utils.gatherchk(y,mode,op.gather);
            
        end % multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Sweepable
            x = matldivide(op,b,mode);
        end % divide
        
    end % Protected methods
    
end % classdef
