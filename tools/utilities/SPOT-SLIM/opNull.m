classdef opNull < opSpot   
%OPNULL  This operator does nothing
%
%   opNull(N) creates the square N-by-N null operator. Without
%   any arguments an operator corresponding to the scalar 1 is
%   created.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        % Constructor
        function op = opNull(n)
            if nargin < 1, n = 1; end
            op = op@opSpot('Null',n,n);
            op.isDirac   = true;
            op.sweepflag = true;
        end

        function A = double(op)
            % double
            A = eye(size(op));
        end

        function result = xtratests(op)
            %XTRATESTS    User defined tests
            %
            % Just a demo here
            x = randn(3);
            A = opNull(3);
            y = A*x;
            result = all(all(x == y));
        end % xtratests

    end % methods - public

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            % Nothing happens here
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,x,mode)
            % and nothing happens here too
        end % divide

    end % methods - protected

end % opNull