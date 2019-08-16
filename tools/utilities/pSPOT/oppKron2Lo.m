classdef oppKron2Lo < oppSpot
    %oppKron2Lo  kronecker tensor product to act on a distributed vector.
    %
    %   oppKron2Lo(A,B, gather)A and B are Spot operators or numeric 
    %   matrices. Optional param gather specifies whether the output vector
    %   should be gathered to local lab.
    %   GATHER = 0 will not gather (default)
    %   GATHER = 1 will gather the results of forwards or adjoint
    %   multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backwards (adjoint) mode.
    %
    %       %% Example: Defining seperable sparsity transforms over 
    %       different axes. Here we define a sparsity transform S that 
    %       performs Wavelet analysis on the first dimension
    %       and a 2D Curvelet analysis on the second & third dimension
    %             dim=[64,32,32];
    %             C = opCurvelet(dim(2),dim(3));
    %             W = opWavelet(dim(1),1);
    %             S = oppKron2Lo(C,W,1);
    %
    %       % Make a random 3d data-array
    %       D = distributed.randn(dim(1),prod(dim(2:end)));
    %
    %       % Check to see if the analysis followed by synthesis returns 
    %       the original signal
    %       norm(D(:)-S'*S*D(:))
    %
    %   note that the second operator is applied to x and then the first
    %   operator to the transpose of the result, and so x should be of
    %   dimemsions [cols(op2),cols(op1)], and vectorized after distribution
    %   so it is distributed along the columns evenly.
    %
    %   *Now oppKron2Lo also supports local x vectors and distributes them
    %   before calculation( this distribution will be faster than using 
    %   matlabs (:) function).
    %
    %   See also: oppKron2Lo.drandn, oppKron2Lo.rrandn, oppKron2Lo.dzeros,
    %   oppKron2Lo.rzeros
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        permutation; %Permutation vector of intergers defining the order to
        %use when the operators (children) of the Kronecker product are
        %applied to a data vector.
        A; % Child operators
        B;
        sizeA;
        sizeB;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppKron2Lo(varargin)
            
            % Setup and extract gather
            if isscalar(varargin{end}) && isnumeric(varargin{end})
                assert(any(varargin{end} == [0 1 2 3]),...
                    'Gather must be 0 1 2 or 3')
                gat = varargin{end};
                varargin(end) = [];
            else
                gat = 0;
            end
            
            % Check for number of operators
            assert(length(varargin) == 2,'Must specify only 2 operators!')
            
            % Standard pSpot Check
            [opList,m,n,cflag,linear] = pSPOT.utils.stdpspotchk(varargin{:});
            opA = opList{1};
            opB = opList{2};
            spmd
                childA = opA;
                childB = opB;
            end
            
            % Construct operator
            op = op@oppSpot('pKron', prod(m), prod(n));
            op.cflag       = cflag;
            op.linear      = linear;
            op.sweepflag   = true;
            op.children    = [];
            op.A           = childA;
            op.B           = childB;
            op.gather      = gat;
            op.permutation = [1 2];
            op.opsn        = n;
            op.opsm        = m;
            op.sizeA       = size(opA);
            op.sizeB       = size(opB);
            
            % Evaluate the best permutation to use when a multiplication is
            % applied
            if ~(prod(m) == 0 || prod(n) == 0)                
                if (m(2)-n(2)) / (m(2)*n(2)) < (m(1)-n(1)) / (m(1)*n(1))
                    op.permutation = [2 1];
                end
            end
            
             % Setting up implicit dimensions of output vector
            % Flipped
            op.ms = {[] []}; op.ns = {[] []};
            if length(opB.ms) > 1
                op.ms{1} = [op.ms{1} opB.ms(:)'];
            else
                op.ms{1} = [op.ms{1} opB.ms{:}];
            end
            if length(opB.ns) > 1
                op.ns{1} = [op.ns{1} opB.ns(:)'];
            else
                op.ns{1} = [op.ns{1} opB.ns{:}];
            end
            
            if length(opA.ms) > 1
                op.ms{2} = [op.ms{2} opA.ms(:)'];
            else
                op.ms{2} = [op.ms{2} opA.ms{:}];
            end
            if length(opA.ns) > 1
                op.ns{2} = [op.ns{2} opA.ns(:)'];
            else
                op.ns{2} = [op.ns{2} opA.ns{:}];
            end
            
        end % Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = headerMod(op,header,mode)
            % Extract explicit size indices
            exsize = header.exsize;
            href   = @spot.data.headerRef; % Used function handles because its shorter
            hasg   = @spot.data.headerAsgn;

            if mode == 1
                
                % Setup variables
                opList = {op.B{1} op.A{1}}; % Last op applied first
                % Number of output dimensions
                n_out_dims = length(spot.utils.uncell(op.ms)) + size(exsize,2) - 1;
                
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
                
                % Replace old first (collapsed) dimensional sizes with operator sizes.
                i = 1;
                x = 1;
                % Extract collapsed first dims from header.
                first_header = href(header,exsize(:,1));
                for u = 1:2
                    % Input header (including collapsed)
                    y            = length(spot.utils.uncell(op.ns{u})) + x - 1;
                    in_header    = href(first_header,x:y);
                    in_header.exsize = [1;y-x+1];
                    % child header
                    child_header = headerMod(opList{u},in_header,mode);
                    % Assignment indices
                    j            = length(spot.utils.uncell(op.ms{u})) + i - 1;
                    % header assignment
                    oldsize      = length(header_out.size);
                    header_out   = hasg(header_out,child_header,i:j);
                    newsize      = length(header_out.size);
                    i            = j + 1 + newsize - oldsize;
                    x            = y + 1;
                end
                exsize_out = 1:length(header_out.size);
                exsize_out = [exsize_out;exsize_out];
                h = header_out;
                h.exsize = exsize_out;
            else
                % Setup variables
                opList = {op.B{1} op.A{1}}; % Last op applied first
                % Number of output dimensions
                n_out_dims = length(spot.utils.uncell(op.ns)) + size(exsize,2) - 1;
                
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
                
                % Replace old first (collapsed) dimensional sizes with operator sizes.
                i = 1;
                x = 1;
                % Extract collapsed first dims from header.
                first_header = href(header,exsize(:,1));
                for u = 1:2
                    % Input header (including collapsed)
                    y            = length(spot.utils.uncell(op.ms{u})) + x - 1;
                    in_header    = href(first_header,x:y);
                    in_header.exsize = [1;y-x+1];
                    % child header
                    child_header = headerMod(opList{u},in_header,mode);
                    % Assignment indices
                    j            = length(spot.utils.uncell(op.ns{u})) + i - 1;
                    % header assignment
                    oldsize      = length(header_out.size);
                    header_out   = hasg(header_out,child_header,i:j);
                    newsize      = length(header_out.size);
                    i            = j + 1 + newsize - oldsize;
                    x            = y + 1;
                end
                exsize_out = 1:length(header_out.size);
                exsize_out = [exsize_out;exsize_out];
                h = header_out;
                h.exsize = exsize_out;
            end
            
        end % headerMod
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            
            childs = {op.A{1} op.B{1}};
            str=['pKron(',char(childs{1})];
            
            % Get operators
            for i=2:length(childs)
                str=strcat(str,[', ',char(childs{i})]);
            end
            str=strcat(str,')');
        end % Char
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double is overloaded to use a vectorized identity matrix
        function y = double(op, enable)
            if nargin < 2 || ~enable
                error(['oppKron2Lo is intended for large applications.' ...
                    ' The explicit representation will likely be very '...
                    ' large, \nuse double(x,1)  to proceed anyway']);
            end
            
            y = double(kron(op.A{1},op.B{1}));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn/rrandn/zeros
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn is overloaded to create a distributed random vector
        function y = drandn(op,varargin)
        %DRANDN Random vector in operator domain
                dims = [op.opsn(2), op.opsn(1)];
            y = distrandnvec( dims );
        end
        function y = rrandn(op,varargin)
            %RRANDN Random vector in operator range
                dims = [op.opsm(2), op.opsm(1)];
            y = distrandnvec( dims );
        end
        function y = dzeros(op,varargin)
            %DZEROS Zero vector in operator domain
                dims = [op.opsn(2), op.opsn(1)];
            y = distzeros( dims );
        end
        function y = rzeros(op,varargin)
            %RZEROS Zero vector in operator range
                dims = [op.opsm(2), op.opsm(1)];
            y = distzeros( dims );
        end
        
    end % Methods
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % The Kronecker product (KP) is applied to the right-hand matrix
            % taking in account the best order to apply the operators A and
            % B.
            
            % Check for op and  x
            assert(isa(op,'oppKron2Lo'),'Left multiplication not taken in account')
            assert( isvector(x) , 'Please use vectorized matrix')
            
            %Operators
            opA = op.A;
            opB = op.B;

            % Transpose checking
            if mode == 2
                op.permutation = fliplr(op.permutation);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%Multiplication%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            if op.permutation(1)==2 %Classic multiplication order
                % Distribute if necessary
                if ~isdistributed(x)
                    if mode == 2
                        sizeA=fliplr(op.sizeA); sizeB=fliplr(op.sizeB);
                    else
                        sizeA=op.sizeA; sizeB=op.sizeB;
                    end
                    rA = sizeA(1); cA = sizeA(2);
                    rB = sizeB(1); cB = sizeB(2);
                    %fprintf('0c:(%d,%d)(%d,%d)\n',rA,cA,rB,cB)
                    y=reshape(x,cB,cA);
                    y=distributed(y);
                else
                    y=x;
                end
                    
                spmd
                    % For transpose case
                    if mode == 2
                        opA = opA';
                        opB = opB';
                    end
                    
                    % Size of the operators
                    [rA,cA] = size(opA);
                    [rB,cB] = size(opB);
                    %fprintf('1c:(%d,%d)(%d,%d)\n',rA,cA,rB,cB)
                    
                    y         = getLocalPart(y);
                    loc_width = prod(size(y))/cB;
                    assert( mod(loc_width,1) == 0, ...
                        'x must be distributed along cols before vec')
                    y         = reshape(y,cB,loc_width); % reshape to local 
                    
                    if ~opB.isDirac
                        y=opB*y;% apply B to local matrices
                    end
                    
                    if ~opA.isDirac
                        y=y.'; % Tranpose
                        % Build y distributed across rows
                        y = codistributed.build(y, codistributor1d...
                            (1,[],[cA,rB]),'noCommunication');
                        y = redistribute(y,codistributor1d(2));%distributed
                        % along cols
                        y = getLocalPart(y);
                        y = opA*y;% apply A to local matrices, then transpose
                        y = y.'; % Now y is distributed across rows again
                        % Rebuild y as distributed across rows
                        y = codistributed.build(y,codistributor1d(1,...
                            [],[rB,rA]),'noCommunication');
                        % Redistribute y across columns
                        y = redistribute(y,codistributor1d(2));
                        y = getLocalPart(y);
                    end
                    
                    % now vectorize y                    
                    % Setup part of columns
                    part = codistributed.zeros(1,numlabs);
                    part(labindex) = numel(y);
                    y              = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        part, [rA*rB,1]),'noCommunication');
                end
            else  %Inverted multiplication order
                % Distribute if necessary
                if ~isdistributed(x)
                    if mode == 2
                        sizeA=fliplr(op.sizeA); sizeB=fliplr(op.sizeB);
                    else
                        sizeA=op.sizeA; sizeB=op.sizeB;
                    end
                    rA = sizeA(1); cA = sizeA(2);
                    rB = sizeB(1); cB = sizeB(2);
                    %fprintf('0i:(%d,%d)(%d,%d)\n',rA,cA,rB,cB)
                    y=reshape(x,cB,cA);
                    y=distributed(y);
                else
                    y=x;
                end
                                        
                spmd
                    % For transpose case
                    if mode == 2
                        opA = opA';
                        opB = opB';
                    end
                    
                    % Size of the operators
                    [rA,cA] = size(opA);
                    [rB,cB] = size(opB);
                    %fprintf('1i:(%d,%d)(%d,%d)\n',rA,cA,rB,cB)

                    y         = getLocalPart(y);
                    loc_width = prod(size(y))/cB;
                    assert( mod(loc_width,1) == 0, ...
                        'x must be distributed along cols before vec')
                    y         = reshape(y,cB,loc_width); % reshape to local 
                    
                    if ~opA.isDirac
                        y = y.'; % transpose since A is applied first
                        % Build y distributed across rows
                        y = codistributed.build(y, codistributor1d...
                            (1,[],[cA,cB]),'noCommunication');
                        % Redistribute y across cols
                        y = redistribute(y,codistributor1d(2));
                    end
                    
                    if ~opA.isDirac
                        y = getLocalPart(y);
                        y = opA*y;% apply A to local matrices, then transpose
                        y = y.';
                        
                        % Rebuild y distributed across rows
                        y = codistributed.build(y,codistributor1d(1,...
                            [],[cB,rA]),'noCommunication');
                        % Redistribute y across cols
                        y = redistribute(y,codistributor1d(2));
                        y = getLocalPart(y);
                    end
                    
                    if ~opB.isDirac
                        y = opB*y;% apply B to local matrices, no need to 
                    end         % transpose
                    
                    % now vectorize y
                    % Setup part of columns
                    part = codistributed.zeros(1,numlabs);
                    part(labindex) = numel(y);
                    y              = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        part, [rA*rB,1]),'noCommunication');
                end
            end
            % if op.gather, y = gather(y); end %#ok<PROP,CPROP>
            if mode == 2 
                if op.gather == 1 || op.gather == 3           
                    y = gather(y);                            
                end                                           
            else % this is the forward case
                if op.gather == 1 || op.gather == 2
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
    end %Protected methods
    
end % Classdef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distrandnvec helper funtion to create a vectorized distributed
% normal random matrix
function y = distrandnvec( sz )
    spmd
        y      = codistributed.randn(sz);
        codist = getCodistributor(y);
        part   = codist.Partition;
        part   = part.*sz(1);
        y      = getLocalPart(y);
        y      = y(:);
        y      = codistributed.build(y, codistributor1d(1,part,...
                 [sz(1)*sz(2),1]));
    end
end

% distzeros helper funtion to create a vectorized distributed
% zeros matrix
function y = distzeros( sz )
    spmd
        y      = codistributed.zeros(sz);
        codist = getCodistributor(y);
        part   = codist.Partition;
        part   = part.*sz(1);
        y      = getLocalPart(y);
        y      = y(:);
        y      = codistributed.build(y, codistributor1d(1,part,...
                 [sz(1)*sz(2),1]));
    end
end

