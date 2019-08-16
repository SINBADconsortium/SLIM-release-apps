classdef (HandleCompatible) oppOrthogonal < oppSpot
%OPPORTHOGONAL   Abstract class for orthogonal operators.
%
%   oppOrthogonal methods:
%     oppOrthogonal - constructor
%     mldivide     - solves Ax=b  via  x=A'b.

%   NOTE: There's no reason to overload @oppSpot/mrdivide because it's simply
%   a wrapper to mldivide.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
    
   methods
      
      function op = oppOrthogonal(type,m,n)
         %opOrthogonal   Constructor for the abstract class. 
         op = op@oppSpot(type,m,n);
      end         
         
      function x = mldivide(op,b)
         %\ (backslash)  x = op\b
         x = op'*b;
      end
            
   end % methods
      
end % classdef
