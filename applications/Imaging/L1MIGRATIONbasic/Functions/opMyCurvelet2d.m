classdef opMyCurvelet2d < opSpot
%OPCURVELET  Two-dimensional curvelet operator.
%
%   opCurvelet(M,N,NBSCALES,NBANGLES,TTYPE) creates a two-dimensional
%   curvelet operator for M by N matrices. The curvelet transform is
%   computed using the Curvelab code.
%
%   The remaining three parameters are optional; NBSCALES gives the
%   number of scales and is set to max(1,ceil(log2(min(M,N)) - 3)) by
%   default, as suggested by Curvelab. NBANGLES gives the number of
%   angles at the second coarsest level which must be a multiple of
%   four with a minimum of 8. By default NBANGLES is set to 16. TTYPE
%   determines the type of transformation and is set to 'WRAP' by
%   default.
%
%   See also CURVELAB.

%   Copyright 2009, Gilles Hennenfent, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

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
       function op = opMyCurvelet2d(m,n,nbscales,nbangles,is_real)

		if nargin < 3, nbscales = ceil(log2(min(m,n)) - 3); end;
		if nargin < 4, nbangles = 16;                       end;
		if nargin < 5, is_real = 1;                         end;

		finest  = 1;

		% Compute length of curvelet coefficient vector
		C = fdct_wrapping_mex(m,n,nbscales,nbangles,finest,randn(m,n));

		hdr{1}{1} = size(C{1}{1});
		cn = prod(hdr{1}{1});  
		for i = 2:nbscales
		    nw = length(C{i});
		    hdr{i}{1} = size(C{i}{1});
		    hdr{i}{2} = size(C{i}{nw/4+1});
		    cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
		end

		parms = {m,n,cn,hdr,finest,nbscales,nbangles,is_real};
		fun = @(x,mode) opMyCurvelet2d_intrnl(parms{:},x,mode);

          % Construct operator
          op = op@opSpot('Curvelet', m*n,cn);
          op.cflag     = ~is_real;
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


function y = opMyCurvelet2d_intrnl(m,n,cn,hdr,ac,nbs,nba,is_real,x,mode)

if mode == 1
    % Synthesis mode  
    x = fdct_v2c(x,hdr,ac,nba);
    if is_real
        x = fdct_wrapping_r2c(x);
        y = real(ifdct_wrapping_mex(m,n,nbs,nba,ac,x));
    else
        y = ifdct_wrapping_mex(m,n,nbs,nba,ac,x);
    end
    y = y(:);
else
    % Analysis mode
    y = fdct_wrapping_mex(m,n,nbs,nba,ac,reshape(x,m,n));
    if is_real
        y = fdct_wrapping_c2r(y);
    end
    y = fdct_c2v(y,cn);
end
end
