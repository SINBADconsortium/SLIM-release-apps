classdef oppDirac < oppOrthogonal   
%OPPDIRAC  Parallel Dirac basis.
%
%   oppDirac(N) creates the square N-by-N identity operator. Without
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
      function op = oppDirac(n)
         if nargin < 1, n = 1; end
         op = op@oppOrthogonal('pDirac',n,n);
         op.isDirac = true;
         op.sweepflag  = true;
      end
      
      function A = double(op)
         A = eye(size(op));
      end
      
      function result = xtratests(op)
      %XTRATESTS    User defined tests
      %
      % Just a demo here
      result = true;
          disp('How thoughtful of you to test oppDirac!!!');
      end
      
   end % methods - public
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - protected
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )
      
      % Multiplication
      function y = multiply(op,x,mode)
         y = x;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Divide
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function x = divide(op,b,mode)
          x = b;
      end % divide
      
   end % methods - protected
   
end % classdef
