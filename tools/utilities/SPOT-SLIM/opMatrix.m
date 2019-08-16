classdef opMatrix < opSpot
%OPMATRIX   Convert a numeric matrix into a Spot operator.
%
%   opMatrix(A,DESCRIPTION) creates an operator that performs
%   matrix-vector multiplication with matrix A. The optional parameter
%   DESCRIPTION can be used to override the default operator name when
%   printed.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        matrix = {}; % Underlying matrix
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opMatrix(A,description)
            % opMatrix  Constructor.
            if nargin < 1
                error('At least one argument must be specified.')
            end
            if nargin > 2
                error('At most two arguments can be specified.')
            end

            % Check if input is a matrix
            if ~(isnumeric(A) || issparse(A))
                error('Input argument must be a matrix.');
            end

            % Check description parameter
            if nargin < 2, description = 'Matrix'; end

            % Create object
            op = op@opSpot(description, size(A,1), size(A,2));
            op.cflag     = ~isreal(A);
            op.sweepflag = true;
            op.matrix    = A;
        end % constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % char
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % char  Create character array from operator.
            if isscalar(op)
                v   = op.matrix;
                str = strtrim(evalc('disp(v)'));
            else
                str = char@opSpot(op);
            end          
        end % function char
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = double(op)
            % double  Convert operator to a double.
            x = op.matrix;
        end % double
        
    end % public Methods


    methods ( Access = protected )

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            % multiply  Multiply operator with a vector.
            if mode == 1
                x = op.matrix * x;
            else
                x = op.matrix' * x;
            end
        end % function multiply

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,x,mode)
            % divide  Solve a linear system with the operator.
            if mode == 1
                x = op.matrix \ x;
            else
                x = op.matrix' \ x;
            end
        end % function divide

    end % private methods
   
end % opMatrix