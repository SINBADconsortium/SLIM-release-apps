classdef opdWindowLast1Halo < oppSpot
%opdWindowLast1Halo tapered windowing for CARP-CG method
%
%   opdWindowLast1Halo(N,L,P,H)
%
%   ARGUMENTS:
%      N = length of the input vector
%      L = length of the last dimension of input vector
%      P = number of processors
%      H = half of the overlap's size
%
%   ATTRIBUTES:
%      xshape = (p,3) vector holding start, size, end indecies
%           of the default distribution of the input vector in every window
%      yshape = (p,3) vector holding start, size, end indecies
%           of the default distribution of the output vector in every window
%
%   Notes:
%       1. This is a parallel/distributed operator
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
    f = 0;
    l = 0;
    p = 0;
    h = 0;
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
       function op = opdWindowLast1Halo(varargin)
      assert(nargin==4,'dWindowLast1Halo: wrong # of arguments')
      n = varargin{1};
      l = varargin{2};
      p = varargin{3};
      h = varargin{4};
      assert(mod(n,l)==0,'Fatal error: L is not a valid last dimension')
      f = n/l;
          [ m xs ys ] = pSPOT.pWindow.funWindowLast1HaloShape( l, p, h );
      op = op@oppSpot('dWindowLast1Halo',f*m,f*l);
      op.f = f;
      op.l = l;
      op.p = p;
      op.h = h;
      op.yshape = ys;
      op.xshape = xs;
       end % function opdWindowLast1Halo
       
       % xtratests
       function result = xtratests(op)
       T = 14;
       x0=distributed.randn(op.f,op.l);
       x0=x0(:);
       x1=op'*op*x0;
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
           y = pSPOT.pWindow.funWindowLast1HaloMakeDist(x,op.l,op.p,op.h);
           else
           y = pSPOT.pWindow.funWindowLast1HaloDropDist(x,op.l,op.p,op.h);
           end
        end % Multipy
      
    end % Methods
        
end % Classdef
