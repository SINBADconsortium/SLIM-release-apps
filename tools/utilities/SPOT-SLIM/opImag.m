classdef opImag < opSpot
%OPIMAG   Complex imaginary part of operator.
%
%   opImag(OP) is the complex imaginary part of operator OP. Note
%   that the resulting operator is real.
%
%   See also opConj, opReal.

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
        function op = opImag(A)

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
            op = op@opSpot('Imag', m, n);
            op.cflag     = false;
            op.linear    = A.linear;
            op.sweepflag = A.sweepflag;
            op.children  = {A};
        end % Constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            op1 = op.children{1};
            str = ['Imag(', char(op1), ')'];
        end % Char

    end % public Methods

    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            opA = op.children{1};
            if mode == 1
                if isreal(x)
                    % Purely real
                    x = imag(applyMultiply(opA,x,mode));
                elseif isreal(1i*x)
                    % Purely imaginary
                    x = imag(applyMultiply(opA,imag(x),mode)) * 1i;
                else
                    % Mixed
                    x = imag(applyMultiply(opA,real(x),mode)) + ...
                        imag(applyMultiply(opA,imag(x),mode)) * 1i;
                end
            else
                if isreal(x)
                    % Purely real
                    x = imag(applyMultiply(opA,x,mode)) * -1;
                elseif isreal(1i*x)
                    % Purely imaginary
                    x = imag(applyMultiply(opA,imag(x),mode)) * 1i * -1;
                else
                    % Mixed
                    x = imag(applyMultiply(opA,real(x),mode)) * -1 + ...
                        imag(applyMultiply(opA,imag(x),mode)) * 1i * -1;
                end
            end
        end % Multiply
       
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

    end % protected Methods
   
end % opImag