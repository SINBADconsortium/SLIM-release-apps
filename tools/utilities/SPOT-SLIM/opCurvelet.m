classdef opCurvelet < opSpot
%OPCURVELET  Two-dimensional curvelet operator.
%
%   opCurvelet(M,N,NBSCALES,NBANGLES,FINEST,TTYPE,IS_REAL) creates a
%   two-dimensional curvelet operator for M by N matrices. The curvelet
%   transform is computed using the Curvelab code.
%
%   The remaining five parameters are optional; NBSCALES gives the number
%   of scales and is set to max(1,ceil(log2(min(M,N)) - 3)) by default, as
%   suggested by Curvelab. NBANGLES gives the number of angles at the
%   second coarsest level which must be a multiple of four with a minimum
%   of 8. By default NBANGLES is set to 16. FINEST sets whether to include
%   the finest scale of coefficients and is set to 0 by default; set this
%   to 1 to include the finest scale, or to 2 to keep the finest scale but
%   set it to zeros. TTYPE determines the type of transformation; either
%   'WRAP' for a wrapping transform or 'ME' for a Mirror-Extended
%   Transform, it's set to 'WRAP' by default.  IS_REAL sets whether the
%   transform is for real data or complex data.
%
%   See also CURVELAB.

%   An example that actually works is opCurvelet(32,32), use this for
%   testing.

%   Nameet Kumar - Oct 2010
%   Copyright 2009, Gilles Hennenfent, Ewout van den Berg and 
%   Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        nbscales;
        nbangles;
        finest;   
        header;          %sizes of coefficient vectors
        nbcoeffs;           %total number of coefficients
        dims;           %size of curvelet
        ttype;           %type of transformation
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opCurvelet(m,n,nbscales,nbangles,finest,ttype,...
                                                                   is_real)

            assert( isscalar(m) && isscalar(n),['Please ensure'...
            ' sizes are scalar values']);
            if nargin < 3, nbscales = max(1,ceil(log2(min(m,n)) - 3)); end;
            if nargin < 4, nbangles = 16;                              end;
            if nargin < 5, finest = 0;                                 end;
            if nargin < 6, ttype = 'WRAP';                             end;
            if nargin < 7, is_real = 1;                                end;
            assert( strcmp(ttype,'WRAP') || strcmp(ttype,'ME'),...
                ['Please ensure ttype is set correctly. Options are'...
                ' "WRAP" for a wrapping transform and "ME" for a'...
                ' mirror-extended transform']);
            assert( isscalar(nbscales) && isscalar(nbangles),...
                'Please ensure nbscales and nbangles are scalar values');
            assert((any(finest == [0 1 2])) && (is_real==0||is_real==1),...
                'Please ensure finest and is_real are appropriate values');
            if finest==0, assert( nbscales>1, ['Please ensure that '...
                    'm and n are large enough for nbscales to be '...
                    'greater than 1 while finest is set to 0']); 
            end

            % Compute length of curvelet coefficient vector
            if strcmp(ttype,'ME')
                C = mefcv2(randn(m,n),m,n,nbscales,nbangles);

                hdr{1}{1} = size(C{1}{1});
                cn = prod(hdr{1}{1});
                for i = 2:nbscales
                    nw = length(C{i});
                    hdr{i}{1} = size(C{i}{1});
                    hdr{i}{2} = size(C{i}{nw/2+1});
                    cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
                end
            else
                [tmphdr, cn] = fdct_sizes_mex(m,n,nbscales,nbangles,...
                                                          logical(finest));
                hdr = cell(1,nbscales);
                hdr{1} = {[tmphdr{1:2}]}; 
                for i = 2:nbscales - (~finest)
                    j = 3 + 5*(i-2);
                    hdr{i}={[tmphdr{j+1:j+2}];[tmphdr{j+3:j+4}];...
                            [tmphdr{j}]};
                end
                if ~finest,  hdr{end} = {[tmphdr{end-1:end}];1}; end;
            end

            % Construct operator
            op = op@opSpot('Curvelet', cn, m*n);
            op.cflag     = ~is_real;
            op.nbscales  = nbscales;
            op.nbangles  = nbangles;
            op.finest    = finest;
            op.header    = hdr;
            op.nbcoeffs  = cn;
            op.dims      = [m,n];
            op.ttype     = ttype;
            op.ns        = {[m n]};
            op.sweepflag = true;
        end % Constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rrandn             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overloaded to produce a vector that really falls in the range of
        % op
        function y = rrandn(op)
            y = op.drandn;
            y = multiply(op,y,1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overloaded to modify metadata correctly
        function h = headerMod(op,header,mode)
            exsize = header.exsize;

            if mode == 1
                h = header; % Copy header
                % Replace old first (collapsed) dimensional sizes with 
                % operator sizes.
                h.size(exsize(1,1):exsize(2,1)) = [];
                h.size   = [op.ms{:} h.size];
                h.origin(exsize(1,1):exsize(2,1)) = [];
                h.origin = [0 h.origin];
                h.delta(exsize(1,1):exsize(2,1)) = [];
                h.delta  = [1 h.delta];
                h.label(exsize(1,1):exsize(2,1)) = [];
                h.label  = ['lcurvelet' h.label];
                h.unit(exsize(1,1):exsize(2,1)) = [];
                h.unit   = ['ucurvelet' h.unit];

            else % mode == 2
                h = header;
                h.size(exsize(1,1):exsize(2,1)) = [];
                h.size   = [op.ns{:} h.size];
                h.origin(exsize(1,1):exsize(2,1)) = [];
                h.origin = [0 0 h.origin];
                h.delta(exsize(1,1):exsize(2,1)) = [];
                h.delta  = [1 1 h.delta];
                h.label(exsize(1,1):exsize(2,1)) = [];
                h.label  = ['l1' 'l2' h.label];
                h.unit(exsize(1,1):exsize(2,1)) = [];
                h.unit   = ['u1' 'u2' h.unit];

            end

            % Re-append correct exsize
            exsize_out = 1:length(h.size);
            exsize_out = [exsize_out;exsize_out];
            h.exsize   = exsize_out;
        end % headerMod
    end % Public Methods
       
 
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            x_n = size(x,2);
            
            if mode == 1
                for u = x_n:-1:1 % Loop over multivectors
                    x_tmp = x(:,u);
                    % Analysis mode
                    if strcmp(op.ttype,'ME')
                        x_tmp = mefcv2(reshape(x_tmp,op.dims(1),...
                                op.dims(2)),op.dims(1),op.dims(2),...
                                op.nbscales,op.nbangles);
                    else
                        x_tmp = fdct_wrapping_mex(op.dims(1),op.dims(2),...
                                op.nbscales,op.nbangles,...
                                logical(op.finest),...
                                reshape(x_tmp,op.dims(1),op.dims(2)));
                        if op.finest == 2, zero_finest_scale; end
                        if ~op.cflag % real transforms have redundancy 
                            x_tmp = fdct_wrapping_c2r(x_tmp);
                        end
                    end
                    y(:,u) = spot.utils.fdct_c2v(x_tmp,op.nbcoeffs);
                end
                x = y;
            else
                for u = x_n:-1:1 % Loop over multivectors
                    x_tmp = x(:,u);
                    % Synthesis mode  
                    if strcmp(op.ttype,'ME')
                        x_tmp = spot.utils.mefdct_v2c(full(x_tmp),...
                            op.header,op.nbangles);
                        x_tmp = meicv2(x_tmp,op.dims(1),op.dims(2),...
                            op.nbscales,op.nbangles);
                    else
                        x_tmp = spot.utils.fdct_v2c(x_tmp,op.header);
                        if op.finest == 2, zero_finest_scale; end
                        if ~op.cflag
                            x_tmp = fdct_wrapping_r2c(x_tmp);
                        end
                        x_tmp = ifdct_wrapping_mex(op.dims(1),...
                                op.dims(2),op.nbscales,op.nbangles,...
                                logical(op.finest),x_tmp);
                    end
                    if ~op.cflag % real transforms don't have complex
                                 % numbers for the inverse transform
                        x_tmp = real(x_tmp);
                    end
                    y(:,u) = x_tmp(:);
                end
                x = y;
            end

            %%% Nested Function
            function zero_finest_scale
                for i = 1:length(x_tmp{end})
                    x_tmp{end}{i} = zeros( size( x_tmp{end}{i} ) );
                end
            end

        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Non-sweepable
            x = lsqrdivide(op,b,mode);
        end % divide

    end % Protected Methods
end % opCurvelet