classdef opMask < opSpot
%OPMASK  Selection mask.
%
%   opMask(N,IDX) creates a diagonal N-by-N operator that has ones
%   only on those locations indicated by IDX.
%
%   opMask(IDX) is the same as opMask(numel(IDX),IDX) when parameter
%   IDX is logical.
%
%   See also opDiag, opRestriction.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        mask = []; % Binary mask
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opMask(varargin)

            if nargin == 1 && islogical(varargin{1})
                idx = varargin{1};
                n   = numel(idx);
            elseif nargin == 2
                n   = varargin{1};
                idx = varargin{2};
            else
                error('Invalid number of parameters specified.')
            end

            idx = full(idx(:));

            if islogical(idx)
                if length(idx) > n
                    error('Index exceeds operator dimensions.');
                end
            elseif spot.utils.isposintmat(idx) || isempty(idx)
                if ~isempty(idx) && (max(idx) > n)
                    error('Index exceeds operator dimensions.');
                end
            else
                error(['Subscript indices must either be real positive '...
                       'integers or logicals.']);
            end

            % Set mask
            mas      = zeros(n,1);
            mas(idx) = 1;

            % Construct operator
            op = op@opSpot('Mask',n,n);
            op.mask      = mas;
            op.sweepflag = true;
        end % Constructor
        
    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function x = multiply(op,x,mode)
            for i=size(x,2):-1:1
                x(:,i) = op.mask.*x(:,i);
            end
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Non-sweepable
            x = lsqrdivide(op,b,mode);
        end % divide
      
    end % Methods
        
end % Classdef
