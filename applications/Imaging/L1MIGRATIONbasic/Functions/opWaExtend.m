classdef opWaExtend < opSpot
%  op = opWaExtend(nz,nx,ttype)
%  OPWAEXTEND is a spot operator which extend the size of image for the
%  Wavetom transform
%  ttype:
%       'pad', padding zeros
%       'extend', extend the model with outer part
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
       function op = opWaExtend(nz,nx,ttype)
		m = 1;
		while 2^m < max(nx,nz)
			m = m + 1;
			n = 2^m;
		end
		
		if nargin < 3, ttype = 'extend';end
		
		fun = @(x,mode) opWaExtend_intrnl(nz,nx,n,ttype,x,mode);

          % Construct operator
          op = op@opSpot('WaveAtomExtend',n^2,nz*nx);
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


function y = opWaExtend_intrnl(nz,nx,n,ttype,x,mode)


if (mode == 1)  % padding zeros
	if strcmp(ttype,'pad')
		x  = reshape(x,nz,nx);
		y  = zeros(n,n);
		sz = round(.5*(n - nz));
		sx = round(.5*(n - nx));
		y(sz+1:sz+nz,sx+1:sx+nx) = x;
		y = y(:);
	elseif strcmp(ttype,'extend')
		sz = round(.5*(n - nz));
		sx = round(.5*(n - nx));
		xextend = opMatrix([ones(sx,1) zeros(sx,nx-1);eye(nx);zeros(n-sx-nx,nx-1) ones(n-sx-nx,1)]);
		zextend = opMatrix([ones(sz,1) zeros(sz,nz-1);eye(nz);zeros(n-sz-nz,nz-1) ones(n-sz-nz,1)]);
		extend  = opKron(xextend,zextend);
		y       = extend * x;
	end
	
elseif (mode == 2)   % remove border
    
	x = reshape(x,n,n);
	sz = round(.5*(n - nz));
	sx = round(.5*(n - nx));
	y = x(sz+1:sz+nz,sx+1:sx+nx);
	y = y(:);

end
end
