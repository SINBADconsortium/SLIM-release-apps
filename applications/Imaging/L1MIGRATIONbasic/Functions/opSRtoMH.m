classdef opSRtoMH < opSpot
% OPSRtoMH  Map source-receiver to midpoint offset.
%
%    OPSRtoMH(A,OPINFO) creates an operator that performs
%    matrix-vector multiplication with matrix A. Optional parameter
%    OPINFO can be used to override the default operator
%    information returned when querying the operator, or provide a 
%    string giving the interpretation of the matrix. When OPINFO is
%    a cell array the first entry is the operator name followed by
%    additional information for that operator type (this option is
%    mostly provided for internal use by other operators).
%
%
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       funHandle = []; % Multiplication function
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opSRtoMH(n)
	
		[XX,YY] = meshgrid(1:n,1:n);
		h       = n+XX-YY;
		m       = XX+YY-1;
		xs      = (m+h+1-n)/2;
		xr      = (m-h+1+n)/2;
		h       = h(:);
		m       = m(:);
		xs      = xs(:);
		xr      = xr(:);
		
		if nargin < 2
		  opinfo = {'opSR2MH', []};
		elseif ischar(opinfo)
		  opinfo = {'opSR2MH', opinfo};
		end

		fun = @(x,mode) opSRtoMH_intrnl(n,x,mode,opinfo,h,m,xs,xr);

          % Construct operator
          op = op@opSpot('SRtoMH',2*n*n,n*n);
          op.cflag     = 1;
          op.funHandle = fun;
       end % Constructor

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = op.funHandle(x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%=======================================================================


function y = opSRtoMH_intrnl(n,x,mode,opinfo,h,m,xs,xr) 


if (mode == 1)  % demigration mode 
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
elseif (mode == 2)   % migration mode
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
end