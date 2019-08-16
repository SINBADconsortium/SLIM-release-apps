classdef opPower < opSpot
%OPPOWER   Raise operator to integer power.
%
%   opPower(OP,P) creates the operator OP^P for integer values of
%   P. When P = 0, the identity matrix is returned. When P < 0 we
%   reformulate the operator as inv(OP^|P|).

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties( SetAccess = private )
        funHandle      % Multiplication function
        exponent       % Exponent of power operator
    end % properties - private

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opPower(A,p)

            if nargin ~= 2
                error('Exactly two operators must be specified.')
            end
            if p ~= round(p)
                error('Second argument to opPower must be an integer.');
            end

            % Input matrices are immediately cast as opMatrix's.
            if isa(A,'numeric'), A = opMatrix(A); end

            % Check that the input operators are valid.
            if ~isa(A,'opSpot')
                error('Input operator is not valid.')
            end

            % Construct operator
            [m, n] = size(A);
            op = op@opSpot('Power', n, m);
            op.cflag    = A.cflag;
            op.linear   = A.linear;
            op.children = {A};
            op.exponent = p;

            % Create function handle
            if p == 0
                fun = @(op,x,mode) x;
            elseif p > 0
                fun = @opPower_intrnl;
            else
                fun = @(op,x,mode) applyMultiply(inv(A^abs(p)),x,mode);
            end
            op.funHandle = fun;
            op.sweepflag = true;
        end % constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            str = [char(op.children{1}),sprintf('^%d',op.exponent)];
        end % Char

    end % methods - public

    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            x = op.funHandle(op,x,mode);
        end % function multiply

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,x,mode)
            % Non-sweepable
            x = lsqrdivide(op,x,mode);
        end % divide

    end % methods - protected

    methods( Access = private )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % opPower_intrnl
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = opPower_intrnl(op,x,mode)
            p = op.exponent;
            A = op.children{1};
            for ind=1:p
                x = applyMultiply(A,x,mode);
            end
        end % function opPower_intrnl

    end % methods - private
    
end % classdef

