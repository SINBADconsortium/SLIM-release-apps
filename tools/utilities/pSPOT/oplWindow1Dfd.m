classdef oplWindow1Dfd < opSpot
%oplWindow1Dfd windowing for finite difference algorithms
%
%   oplWindow1Dfd(N,P,H)
%
%   ARGUMENTS:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%
%   ATTRIBUTES:
%      oshape = (p,3) vector holding start, size, end indecies
%           of the default distribution of the input vector in every window
%      yshape = (p,3) vector holding start, size, end indecies
%           of the default distribution of the output vector in every window
%      xshape = (p,3) vector holding start, size, end indecies
%           of the input vector that will end up in every ouput window
%
%   Notes:
%       1. This is not a parallel/distributed operator
%       2. This operator does not pass dottest
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        p = 0;
    h = 0;
    oshape = 0;
    yshape = 0;
    xshape = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = oplWindow1Dfd(varargin)
      assert(nargin==3,'lWindow1Dfd: wrong # of arguments')
      n = varargin{1};
      p = varargin{2};
      h = varargin{3};
          [ m os ys xs ] = pSPOT.pWindow.funWindow1DShape( n, p, h );
      op = op@opSpot('lWindow1Dfd',m,n);
      op.p = p;
      op.h = h;
      op.oshape = os;
      op.yshape = ys;
      op.xshape = xs;
       end % function oplWindow1Dfd
       
       % xtratests
       function result = xtratests(op)
           T = 14;
       x0=rand(op.n,1);
       y=op*x0;
       x1=op'*y;
       check=norm(x1-x0);
       if check < op.n*10^-T
               result = true;
       else
               result = false;
       end
       end % xtratests

       % utest ( skip dottest )
       function output = utest(op,k,verbose)
           try
               addpath(fullfile(spot.path,'tests','xunit'))
           catch ME
               error('Can''t find xunit toolbox.')
           end
           if nargin < 3, verbose = 0; end
           if nargin < 2, k = 5; end
           assertTrue(op.xtratests,k);
           output = 'PASSED!';
       end % utest

    end % Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if (mode == 1)
          [ A ] = pSPOT.pWindow.funWindow1DfdFor(op.n,op.p,op.h);
              y = A * x;
           else
          [ B ] = pSPOT.pWindow.funWindow1DfdBck(op.n,op.p,op.h);
              y = B * x;
           end
        end % Multipy
      
    end % Methods
        
end % Classdef
