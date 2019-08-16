classdef oppKron < oppSpot
    %OPPKRON    Kronecker tensor product to act on a data container.
    %   we need an oppKron that can be robust about the distribution
    %   dimension.
    %   oppKron(OP1,OP2,...,OPN) Where OPi are Spot operators or
    %   numeric matrices. Note that the last operator is applied to x then
    %   the rest of the operators to the tranpose of the result, ans so x
    %   should be of dimensions [cols(OP1),cols(OP2),...,cols(OPN)], and
    %   vectorized after distribution.
    %
    %   Optional parameter gather specifies whether the output vector
    %   should be gathered to the local lab.
    %   GATHER = 0 will leave them distributed.
    %   GATHER = 1 will gather the results of forwards or adjoint multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backward (adjoint) mode.
    %
    %   See also: oppKron2Lo
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppKron(varargin)
            % Check parallel pool
            if parpool_size() == 0
                error('Parallel pool is not on');
            end
            
            % Settin' up the variables
            gather = 0;            
            % Extract gather
            if isscalar(varargin{end}) && ~isa(varargin{end},'opSpot')
                gather = varargin{end};
                varargin(end) = [];
            end
            
            % Standard checking and setup sizes
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            m = prod(m);
            n = prod(n);
            
            % Construct operator
            op = op@oppSpot('pKron', m, n);
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
            % Initialize
            str = 'pKron(';
            if ~isnumeric( op.children{1} )
                for child = op.children
                    str = [str,char(child{1}),', '];
                end
            else
                [m,n] = cellfun(@size, op.children);
                for i=1:length(op.children)
                    str = [str,'Matrix(',int2str(m(i)),',', ...
                        int2str(n(i)),'), '];
                end
            end
            str = [str(1:end-2), ')'];
        end % Display
        
    end % methods
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % Check for dataconness of x
            assert(isa(x,'dataContainer'),...
                'X must be a data container')
            assert(x.isdist,'X must be distributed');
            
            % Remove implicit vectorization
            x     = univec(x);
            
            % Setup variables
            ops   = op.children;
            
            % Reversing order of children for intuitive indexing
            ops   = fliplr(ops);
            nops  = length(ops);
            
            % Setup spmd variables
            data  = x.data;
            perm  = circshift(1:nops,[0 -1]);
            ddims = x.codist.Dimension;
            
            % Quickfix for bordercases
            temp  = ones(1,nops);
            temp(1:length(x.dims)) = x.dims;
            gsize = temp;
            
            spmd
                % Setup local parts
                dloc = getLocalPart(data);
                
                % Loop through the children
                for i = 1:nops
                    
                    if i == ddims
                        % Distributed dimension on first dimension
                        % Rebuild and redistribute to second dimension
                        part = codistributed.zeros(1,numlabs);
                        if ~isempty(dloc)
                            part(labindex) = size(dloc,1);
                        end%,labBarrier, fprintf('A: '), size(dloc), part, gsize
                        pcod = codistributor1d(1,part,gsize);
                        dloc = codistributed.build(dloc,pcod,'noCommunication');
                        dloc = redistribute(dloc,codistributor1d(2));
                        dloc = getLocalPart(dloc);
                    end
                    
                    % Multiply
                    if mode == 1
                        OP = ops{i};
                    else
                        OP = ops{i}';
                    end
                    dloc = spot.utils.nDimsMultiply(OP,dloc);
                    
                    % Update gsize
                    gsize(1) = size(OP,1);
                    
                    if i == ddims
                       % Rebuild and redistribute back to first dimension
                       part = codistributed.zeros(1,numlabs);
                       if ~isempty(dloc)
                            part(labindex) = size(dloc,2);
                       end%,labBarrier, fprintf('B: '), size(dloc), part, gsize
                       pcod = codistributor1d(2,part,gsize);
                       dloc = codistributed.build(dloc,pcod,'noCommunication');
                       dloc = redistribute(dloc,codistributor1d(1));
                       dloc = getLocalPart(dloc);
                    end                    
                    
                    % Permute & update gsize
                    dloc  = permute(dloc,perm);
                    gsize = circshift(gsize,[0 -1]);
                end % loop
                
                % Rebuild data
                part = codistributed.zeros(1,numlabs);
                if ~isempty(dloc)
                    part(labindex) = size(dloc,ddims);
                end%,labBarrier, fprintf('C: '), size(dloc), part, gsize
                cod  = codistributor1d(ddims,part,gsize);
                data = codistributed.build(dloc,cod,'noCommunication');
            end % spmd
            
            % Setup the variables
            y = dataContainer(data);
            y.data    = data;
            y.dims    = gsize{1};
            y.codist  = cod{1};
            y.history = x.history;
            setHistory(y);
            y = ivec(y);

        end % Multiply
        
    end % Methods
    
end % classdef










