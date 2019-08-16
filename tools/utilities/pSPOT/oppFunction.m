classdef oppFunction < oppSpot
%OPPFUNCTION   Wrapper for functions.
%
%   oppFunction(M,N,FUN,CFLAG,LINFLAG) creates a wrapper for function
%   FUN, which corresponds to an M-by-N operator. The FUN parameter
%   can be one of two types:
%
%   1) A handle to a function of the form FUN(X,MODE), where the
%      operator is applied to X when MODE = 1, and the transpose is
%      applied when MODE = 2;
%   2) A cell array of two function handles: {FUN,FUN_TRANSPOSE},
%      each of which requires only one parameter, X.
%
%   Optional arguments CFLAG and LINFLAG indicate whether the
%   function implements a complex or real operator and whether it
%   is linear or not. The default values are CFLAG=0, LINFLAG=1.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( SetAccess = private )
       funHandle  % Function handles
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % Constructor
        function op = oppFunction(m,n,funhandle,cflag,linflag)
           if nargin < 3
              error('opFunction requires at least three parameters.');
           end
           if nargin < 4 || isempty(cflag)
              cflag = 0;
           end
           if nargin < 5 || isempty(linflag)
              linflag = 1;
           end
           if ~spot.utils.isposintscalar(m) || ~spot.utils.isposintscalar(n)
              error('Dimensions of operator must be positive integers.');
           end
           
           if iscell(funhandle) && length(funhandle) == 2
              if ~isa(funhandle{1},'function_handle') || ...
                 ~isa(funhandle{2},'function_handle')
                 error('Invalid function handle specified.');
              end
              fun = @(x,mode) oppFunction_intrnl(funhandle,x,mode);
              
           elseif isa(funhandle,'function_handle')
              fun = @(x,mode) funhandle(x,mode);
              
           else
              error('Invalid function handle specified.');
              
           end
           
           % Construct operator
           op = op@oppSpot('pFunction',m,n);
           op.cflag  = cflag;
           op.linear = linflag;
           op.funHandle = fun;
          op.sweepflag  = true;
        end % Constructor
        
    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
            n = size(x,2);
            if isscalar(op)
                % special case: allocate result size of x
                m = size(x,1);
            elseif mode == 1
                m = size(op,1);
            else
                m = size(op,2);
            end
            if isdistributed(x)
                spmd, c = class(getLocalPart(x));end
                y = distributed.zeros(m,n,c{1});
            else
                y = zeros(m,n,class(x));
            end
            
            for u = 1:size(x,2)
                y(:,u) = op.funHandle(x(:,u),mode);
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

%======================================================================

function y = oppFunction_intrnl(funhandle,x,mode)
if mode == 1
   fun = funhandle{1};
else
   fun = funhandle{2};
end

% Evaluate the function
y = fun(x);
end
