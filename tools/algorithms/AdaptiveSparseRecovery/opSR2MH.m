function op = opSR2MH(n)
% OPMATRIX  Map source-receiver to midpoint offset.
%
%    OPMATRIX(A,OPINFO) creates an operator that performs
%    matrix-vector multiplication with matrix A. Optional parameter
%    OPINFO can be used to override the default operator
%    information returned when querying the operator, or provide a 
%    string giving the interpretation of the matrix. When OPINFO is
%    a cell array the first entry is the operator name followed by
%    additional information for that operator type (this option is
%    mostly provided for internal use by other operators). 

%   Copyright 2007, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opMatrix.m 593 2007-09-20 06:04:40Z ewout78 $

% Define the offset operator
[XX,YY] = meshgrid(1:n,1:n);
h       = n+XX-YY;
m       = XX+YY-1;
xs      = (m+h+1-n)/2;
xr      = (m-h+1+n)/2;
h       = h(:);
m       = m(:);
xs      = xs(:);
xr      = xr(:);

% do the map from source-receiver to midpoint-offest

if nargin < 2
  opinfo = {'opSR2MH', []};
elseif ischar(opinfo)
  opinfo = {'opSR2MH', opinfo};
end

op = @(x,mode) opSR2MH_intrnl(n,x,mode,opinfo,h,m,xs,xr);

function y = opSR2MH_intrnl(n,x,mode,opinfo,h,m,xs,xr)
if mode == 0
   y = {2*n*n,n*n,[0,1,0,1],opinfo};
elseif mode == 1
  x = reshape(x,n,n);
  buf      = zeros(2*n,2*n);
  for i=1:n*n
    buf(m(i),h(i)) =x(xs(i),xr(i));
  end
  
  y = zeros(n,2*n);
  y(:,1)=buf(2:2:end,1);
  y(:,end)=buf(1:2:end,end);
  for icol=1:n-1
    y(:,icol*2:icol*2+1)=[buf(1:2:end,icol*2) buf(2:2:end,icol*2+1)];
  end
  y = y(:);
else
  x = reshape(x,n,2*n);
  buf = zeros(2*n,2*n);
  for icol=1:n-1
     buf(1:2:end,icol*2)   =  x(:,icol*2);
     buf(2:2:end,icol*2+1) =  x(:,icol*2+1);
  end
  buf(2:2:end,1)=x(:,1);
  buf(1:2:end,end)=x(:,end);
  y      = zeros(n,n);
  for i=1:n*n
    y(xs(i),xr(i)) = buf( m(i),h(i));
  end
  y = y(:);
end
