classdef oppNumBlockDiag < oppSpot
%OPPNUMBLOCKDIAG   Operator-diagonal operator in parallel sans overlap
%   Supports distributed x and distributed 3D Matrices
%
%   B = oppNumBlockDiag([WEIGHT],A,GATHER) where A is a 3D numerical
%   matrix. This will slice A along the 3rd dimension and use the 2D
%   slices as the blocks for the block-diagonal operator B.
%   Note: 2D x vector not supported for numeric 3D
%   GATHER specifies whether to gather the results to a local array
%   or leave them distributed, default is 0.
%   GATHER = 0 will leave them distributed.
%   GATHER = 1 will gather the results of forwards or adjoint 
%              multiplication.
%   GATHER = 2 will gather only in forward mode.
%   GATHER = 3 will gather only in backward (adjoint) mode.
%
%   B = oppNumBlockDiag([WEIGHT],A,ncols_x,GATHER) where A is a 3D numerical
%   matrix. This will slice A along the 3rd dimension and use the 2D
%   slices Ak to build Kronecker operator kron(opDirac(ncols_x), opMatrix(Ak)).
%   These will then be used as the blocks for the block-diagonal operator B.
%   In this mode the flag GATHER is REQUIRED to distinguish this from the
%   previous case.
%
%
%   See also oppBlockDiag, oppDictionary, oppStack

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

% Tim's change Feb 1, 2011: Now works correctly when number of workers 
% exceed number of operators!! Damn you thomas T_T I spent 2 hours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    locm; % Local sizes for use when multiplying
    locn;
    ncols; % number of implied RHS when handling x that is implicitly 
           % multidimensional, default to 1

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function op = oppNumBlockDiag(varargin)

        % Check parallel pool
        assert(parpool_size() > 0, 'Parallel pool is not open');

        % Setting up the variables
        ncols_x = 1;
        gather  = 0;

        % Arguments Extraction
        % Extract gather & ncols
        if isscalar(varargin{end}) % Gather
            if isscalar(varargin{end-1}) % ncols
                ncols_x         = varargin{end-1};
                varargin(end-1) = [];
            end
            gather        = varargin{end};
            varargin(end) = [];
        end

        % Extract weights
        nargs = length(varargin);
        if isdistributed(varargin{1}) || ~isnumeric(varargin{1}) ||...
           ndims(varargin{1})==3   % No weights
            weights = ones(nargs,1);
        else
            % weights
            weights = varargin{1};
            if isempty(weights), weights = 1; end;
            varargin(1) = []; % Remove
        end

        % Check weights
        nSlices = size(varargin{1},3);
        Cube    = double(varargin{1});
        if isscalar(weights)
            weights = weights.*ones(size(varargin{1},3),1); % default
        else
            if length(weights(:)) == nSlices
                weights = ones(nSlices,1).*weights(:);
            else
                error('WEIGHTS size mismatch');
            end
        end

        % Operators Extraction and Processing
        % Check for distribution of the 3D Matrix
        assert(isdistributed(Cube),'Cube must be distributed!');
        clear varargin;

        % Check complexity 
        [m,n,nSlices]     = size(Cube);
        m = m*nSlices;  n = n*nSlices;
        cflag  = ~isreal(Cube) || ~all(isreal(weights));
        linear = 1;
        % ncols_x for multiple RHS
        m = m*ncols_x;    n = n*ncols_x;


        % Post-processing and Construction of Operator
        % Convert distributed attributes to non-distributed scalars
        spmd
            % Check distribution dimension
            num_cod    = getCodistributor(Cube);            
            glo_ind    = globalIndices(num_cod,3); % Glo inds for weights
            
            % Weights processing
            loc_wgts   = weights(glo_ind);
            % Gather weights
            wgts_size  = codistributed.zeros(1,numlabs);
            wgts_size(labindex) = length(loc_wgts);
            wgts_cod   = codistributor1d(1,wgts_size,[sum(wgts_size) 1]);
            loc_wgts   = codistributed.build(loc_wgts(:),wgts_cod);

            % Extract local m and n for future use
            loc_childs = getLocalPart(Cube);
            [mm,nn,ss] = size(loc_childs);
            loc_m      = mm*ss;
            loc_n      = nn*ss;
        end
        clear Cube;
        num_cod = num_cod{1};
        assert(num_cod.Dimension == 3,...
            '3D Matrix is not distributed correctly!');
        
        % Construct operator
        op = op@oppSpot('pnumBlockDiag', m, n);
        op.locm      = [loc_m{:}]; % Local sums for operator m
        op.locn      = [loc_n{:}]; % Local sums for operator n
        op.cflag     = cflag;
        op.linear    = linear;
        op.children  = loc_childs; % Composite
        op.weights   = loc_wgts; % Distributed
        op.sweepflag = false;
        op.gather    = gather;
        op.ncols     = ncols_x;
        op.opsn      = (n/nSlices)*ones(1,nSlices);
        op.opsm      = (m/nSlices)*ones(1,nSlices);

    end %Constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function str = char(op)
        % Initialize
        str = 'pNumBlockDiag(';
        
        for i=1:length(op.opsn)
            str = strcat(str,'Matrix(',int2str(op.opsm(i)),',',...
                  int2str(op.opsn(i)),'), ');
        end

        % Concatenating composite string
        str = [str(1:end-1) ')'];
    end % Display

end % Methods


methods ( Access = protected )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multiply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y = multiply(op,x,mode)

        % Setting up class variables
        loc_wgts   = op.weights; % is composite
        loc_childs = op.children; % also composite
        n_cols     = op.ncols;

        % Matrix x preprocessing
        assert(size(x,2) == 1, 'Please vectorize your x');
        
        % Forward mode (mode = 1)
        if mode == 1
            % Distributing x into the correct chunks
            x_size = [op.n 1]; % Setting up the global sizes
            y_size = [op.m 1];
            x_part = op.locn;
            y_part = op.locm;
        else % mode == 2
            x_size = [op.m 1];
            y_size = [op.n 1];
            x_part = op.locm;
            y_part = op.locn;
        end % end of mode 2

        % take into account implied multiple RHS in multidimensional x
        x_part = x_part .* n_cols;
        y_part = y_part .* n_cols;

        % Check for distribution of x and redistribute if necessary
        if isdistributed(x)
            warnx = {0};
            spmd
                x_ocod = getCodistributor(x);
                x_ncod = codistributor1d(1,x_part,x_size);
                if x_ocod.Partition ~= x_ncod.Partition
                    warnx = 1;
                    x     = redistribute(x,x_ncod);
                end
            end
            if warnx{1}
                error('x is not distributed correctly!');
            end
            
        else    % Distribute x
            error('x is not distributed!');
        end

        % Multiplication starts
        spmd
            % Multiply using opBlockDiag locally
            loc_x  = getLocalPart(x);
            loc_wgts = getLocalPart(loc_wgts);

            if ~isempty(loc_childs)
                num_ops = size(loc_childs,3); % number of blocks local to 
                                                % this machine
                % reshape x if we have multiple RHS
                % assertion: this block only executes for 3D matrix input. 
                % Therefore, all the sizes (m and n) of all
                % child operators must be the same at this point.
                loc_x = reshape(loc_x, [], n_cols, num_ops);
                if mode == 1
                    y = zeros(size(loc_childs,1), n_cols, num_ops, class(loc_x));
                    for k = 1:num_ops;
                        A        = loc_childs(:,:,k);
                        y(:,:,k) = loc_wgts(k) .* (A * loc_x(:,:,k));
                    end
                else
                    y = zeros(size(loc_childs,2), n_cols, num_ops, class(loc_x));
                    for k = 1:num_ops;
                        A        = loc_childs(:,:,k);
                        y(:,:,k) = conj(loc_wgts(k)).*(A' * loc_x(:,:,k));
                    end
                end

                y = y(:);

            else
                % local part of distributed vector for nodes that do not 
                % contain data is of size (0,1) for some reason
                y = zeros(0,1, class(loc_x));
            end

            % Check for sparsity
            aresparse           = codistributed.zeros(1,numlabs);
            aresparse(labindex) = issparse(y);
            % labBarrier;
            if any(aresparse), y = sparse(y); end;

            y_cod = codistributor1d(1,y_part,y_size);
            y     = codistributed.build(y,y_cod,'noCommunication');

        end % spmd
        
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
            % Non-Sweepable
            x = lsqrdivide(op,b,mode);
        end % divide
end % Protected Methods    
end % Classdef
