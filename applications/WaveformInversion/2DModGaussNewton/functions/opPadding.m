classdef opPadding < opSpot
%  op = opPadding(nz,nx,nzb,nxb)
%  OPPADDING is a spot operator which can padding zeros to a matrix
%  Input:
% 		[nz,nx]: size of the matrix before padding zeros.
% 		[nzb,nxb]: size of the matrix after padding.
%  
%  Forward operator means padding zeros, while adjoint of it means unpadding
%
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


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
       function op = opPadding(nz,nx,nzb,nxb)
	
		% [nz,nx]  = size(opts.Initial_model);
		vb                 = nzb;%round(nz * opts.Vertical_border_rate);
		lb                 = nxb;%round(nx * opts.Lateral_border_rate);
		nzt                = nz + 2 * vb;
		nxt                = nx + 2 * lb;
		fun = @(x,mode) opPadding_intrnl(nz,nx,vb,lb,nzt,nxt,x,mode);

          % Construct operator
          op = op@opSpot('Padding',nzt*nxt,nz*nx);
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


function y = opPadding_intrnl(nz,nx,vb,lb,nzt,nxt,x,mode)


if (mode == 1)  % padding zeros
	if numel(x) == nz*nx;
		x = reshape(x,nz,nx);
		x = [zeros(vb,2*lb + nx); zeros(nz,lb),x,zeros(nz,lb); zeros(vb,2*lb + nx)];
		y = x(:);
	else
		error(['number of elements in the model should be ',num2str(nz*nx)])
	end
elseif (mode == 2)   % remove border
    if numel(x) == nxt*nzt
		x = reshape(x,nzt,nxt);
		x = x(vb+1:vb+nz,lb+1:lb+nx);
		y = x(:);
	else
		error(['number of elements in the model should be ',num2str(nzt*nxt)])
	end

end
end
