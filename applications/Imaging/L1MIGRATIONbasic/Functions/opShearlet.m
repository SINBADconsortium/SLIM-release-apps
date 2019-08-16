classdef opShearlet < opSpot
% OPHSHEARLET  Shearlet operator
%
%    OPSHEARLET(qmf1,qmf2,scale,ndir,h_v) creates a Shearlet
%    operator of  Extended Discrete Shearlet Transform associated
%    with only one type of shear matrix [1 1; 0 1] (or [1 0; 1 1])
%    for M by N matrices. 
%       ( compute shearlet coefficients to produce 
%         directional components so that each of those 
%         components contains shearlet coefficients 
%         associated with the corresponding shearing )
%
%    The shearlet transformation is computed using the shearlab.
%
% Input:
%
%   For j = 1,...,L (where L = size of input vector 'scale'),
%   qmf1: 1D Quadrature mirror filter associated with scaling 2^(L-j) 
%   qmf2: 1D Quadrature mirror filter associated with scaling 2^(L-s(j)) 
%   ( see MakeONFilter.m in WaveLab )
%
%   scale : row vector consisting of sizes of shearlets across scales 
%           scale =  [s(1),...,s(L)] (see sampled_DST.m)
%
%   ndir:   number of directions = 2^(ndir+1)+1;
%   x_y:    
%       1) if x_y = 0 --->
%          2^(ndir+1)+1 directional components assoicated with shear matrices 
%          [1 k/2^(ndir); 0 1] for k = -2^(ndir),...,2^(ndir)
%       2) if x_y = 1 --->
%          2^(ndir+1)+1 directional components assoicated with shear matrices 
%          [1 0; k/2^(ndir) 1] for k = -2^(ndir),...,2^(ndir)
%
%   Author: Xiang Li, SLIM, EOS, UBC


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
       function op = opShearlet(m,n,nbscales,nbangles,is_real,ttype)

		if (nargin < 3), qmf1  = MakeONFilter('Symmlet',4);   end;
		if (nargin < 4), qmf2  = MakeONFilter('Symmlet',4);   end;
		if (nargin < 5), scale  = [3 3 3 4 4];  end;
		if (nargin < 6), ndir   =     0;        end;
		if (nargin < 7), h_v    =     0;        end;


		fun = @(x,mode) opShearlet_intrnl(m,n,qmf1,qmf2,scale,ndir,h_v,x,mode);

          % Construct operator
          op = op@opSpot('Shearlet',n*m,n*m*(2^(ndir+1)+1));
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


function y = opShearlet_intrnl(m,n,qmf1,qmf2,scale,ndir,h_v,x,mode)


if mode == 1
	x     = reshape(x,n,m,(2^(ndir+1)+1));
	if   isreal(x)
		y     = sqrt(1/(2^(ndir+1)+1)).*half_adj_sampled_DST1(x,qmf2,qmf2,scale,ndir,h_v);
	else
		xreal = real(x);
		ximag = imag(x);
		yreal = sqrt(1/(2^(ndir+1)+1)).*half_adj_sampled_DST1(xreal,qmf2,qmf2,scale,ndir,h_v);
		yimag = sqrt(1/(2^(ndir+1)+1)).*half_adj_sampled_DST1(ximag,qmf2,qmf2,scale,ndir,h_v);
		y     = yreal + 1i.*yimag;
	end
	y         = y(:);
	
elseif mode == 2
	x     = reshape(x,n,m);
	if isreal(x)
		y     = sqrt(1/(2^(ndir+1)+1)).*half_sampled_DST1(x,qmf2,qmf2,scale,ndir,h_v);
	else 
		xreal = real(x);ximag = imag(x);
		yreal = sqrt(1/(2^(ndir+1)+1)).*half_sampled_DST1(xreal,qmf2,qmf2,scale,ndir,h_v);
		yimag = sqrt(1/(2^(ndir+1)+1)).*half_sampled_DST1(ximag,qmf2,qmf2,scale,ndir,h_v);
		y     = yreal + 1i.*yimag;
	end
	y         = y(:);
end
end
