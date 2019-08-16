classdef oplWindow1Dtpr < opSpot
%oplWindow1Dtpr tapered windowing for partition-of-unity algorithms
%
%   oplWindow1Dtpr(N,P,H)
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
       function op = oplWindow1Dtpr(varargin)
      assert(nargin==3,'lWindow1Dtpr: wrong # of arguments')
      n = varargin{1};
      p = varargin{2};
      h = varargin{3};
          [ m os ys xs ] = pSPOT.pWindow.funWindow1DShape( n, p, h );
      op = op@opSpot('lWindow1Dtpr',m,n);
      op.p = p;
      op.h = h;
      op.oshape = os;
      op.yshape = ys;
      op.xshape = xs;
       end % function oplWindow1Dtpr
       
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

    end % Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if (mode == 1)
          [ A ] = pSPOT.pWindow.funWindow1DtprFor(op.n,op.p,op.h);
              y = A * x;
           else
          [ B ] = pSPOT.pWindow.funWindow1DtprBck(op.n,op.p,op.h);
              y = B * x;
           end
        end % Multipy
      
    end % Methods
        
end % Classdef
