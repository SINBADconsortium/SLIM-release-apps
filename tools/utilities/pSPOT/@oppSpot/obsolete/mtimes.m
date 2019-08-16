function y = mtimes(A,B)
%*   Product of two operators.
%
%   A*B  returns an operator that is the product of two operators.
%
%   See also opFoG.

%   Nameet Kumar - Oct 2010
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.
   
%   http://www.cs.ubc.ca/labs/scl/spot

% Note that either A or B must be belong to the oppSpot class because this
% function gets called for both M*C, C*M, where C is the class and M is a
% matrix or vector. This gives the following options, with s for scalar and
% C for any instance of an oppSpot:
%
% 1) M*C, implemented as (C'*M')'
% 2) C*M
% 3) s*C
% 4) C*s
% 5) C*C, either of which can be a foreign class (including opSpot)

% dataContainer preprocessing
if isa(B,'dataContainer')
    y = mtimes(B, A,'swap');
else
    
% if nargin == 3 && strcmp(swap,'swap')
%     tmp = A;    A = B;      B = tmp;    clear tmp;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 1: M*C
% Mode 3: s*C - Here we also handle the special case where C is 1-by-M.
%               If so, then we recast this as (C'*s)', which results in
%               a call to the "usual" matrix-vector product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isa(A,'opSpot')
    if isscalar(A) && (B.m ~= 1)
       % s*C (mode 3)
       y = opFoG(A,B);
    else
       % M*C (mode 1)
       y = (B' * A')';
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 2: C*M
% Mode 4: C*s - Here we also handle the special case where C is N-by-1.
%               If so, then we recast this as (C'*s)', which results in
%               a call to the "usual" matrix-vector product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isa(B,'opSpot')
   if isscalar(B)
      if A.n ~= 1
         % C*s (mode 4)
         y = opFoG(A,B);
      else
         y = A.applyMultiply(B,1);  % A is a column "vector".
      end
   else
      p = size(B,1);
   
      % Raise an error when the matrices do not commute. We make an
      % exception for 1-by-1 operators.
      if A.n ~= p
         if A.m == 1 && A.n == 1op.counter.plus1(mode);

                % For border case: empty x
                if isempty(x)
                    if mode == 1
                        y = zeros(op.m,0);
                    else
                        y = zeros(op.n,0);
                    end
                    return
                end

                if op.sweepflag
                    y = op.multiply(x,mode);
                else
                    q = size(x,2);

                    % Preallocate y
                    if q > 1
                        if isscalar(op)
                            % special case: allocate result size of x
                            y = zeros(size(x));
                        elseif mode==1
                            y = zeros(op.m,q);
                        else
                            y = zeros(op.n,q);
                        end
                    end

                    if q == 1
                        y = op.multiply(x,mode);
                    else
                        for i=1:q
                            y(:,i) = op.multiply(x(:,i),mode);
                        end
                    end
                end
            % relax
         else
            error(...
               'Matrix dimensions must agree when multiplying by %s.',...
               char(A));
         end
      end
   
      % Perform operator*matrix
      if isempty(A)
         y = zeros(A.m,0);
      else
         y = A.applyMultiply(B,1);
      end
      
   end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both args are Spot/pSpot ops. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    y = opFoG(A,B);
end
end % catch
