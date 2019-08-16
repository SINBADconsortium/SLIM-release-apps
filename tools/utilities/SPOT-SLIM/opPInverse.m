classdef opPInverse < opSpot
%OPPINVERSE   Pseudo inverse of operator.
%
%   Apinv = opPInverse(A) creates the pseudo inverse of a M-by-N
%   operator A. The product Apinv*b is then equivalent to A\b.
%
%   See also @opSpot/mldivide.

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
        function op = opPInverse(A)

            if nargin ~= 1
                error('Exactly one operator must be specified.')
            end

            % Input matrices are immediately cast as opMatrix's.
            if isa(A,'numeric'), A = opMatrix(A); end

            % Check that the input operators are valid.
            if ~isa(A,'opSpot')
                error('Input operator is not valid.')
            end

            % Construct operator
            [m, n] = size(A);
            op = op@opSpot('PInverse', n, m);
            op.cflag     = A.cflag;
            op.linear    = A.linear;
            op.sweepflag = A.sweepflag;
            op.children  = {A};
            op.ms        = A.ns;
            op.ns        = A.ms;
        end % constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            str = ['pinv(', char(op.children{1}) ,')'];
        end % function char

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PInv
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function opOut = pinv(op)
            opOut = op.children{1};
        end % function pinv

    end % methods - public


    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            opA = op.children{1};
            if mode == 1
                A = opA;
            else
                A = opA';
            end
            x = A\x;
        end % function multiply

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

    end % methods - protected
   
end % opPInverse