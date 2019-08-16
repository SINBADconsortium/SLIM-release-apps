classdef opScale < opSpot
%OPCURVELET  Two-dimensional curvelet operator.
%
%   opCoreScale(M,N,NBSCALES,NBANGLES,FINEST,SUBSCALE,TTYPE) is mask for
%   the Curvelet coefficients which keep the central scales. it's a bit like
%   a low pass filter. 
% 
%   Input arragments refer to spot operator "opCurvelet". Parameters should
%   be exactly the same as the parameters used in Creating Curvelet transform
%   operator.
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
       nbscales;
       nbangles;
       finest; 
  	   subscale;
       header;          %sizes of coefficient vectors
       nbcoeffs;           %total number of coefficients
       dims;           %size of curvelet
       ttype;           %type of transformation
       
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opScale(m,n,nbscales,nbangles,finest,subscale,ttype)

          assert( isscalar(m) && isscalar(n),['Please ensure'...
            ' sizes are scalar values']);
          if nargin < 3, nbscales = max(1,ceil(log2(min(m,n)) - 3)); end;
          if nargin < 4, nbangles = 16;                              end;
          if nargin < 5, finest   = 0;                               end;
		  if nargin < 6, subscale = 1;                               end;
          if nargin < 7, ttype    = 'WRAP';                          end;
		  
          assert( strcmp(ttype,'WRAP') || strcmp(ttype,'ME'), ['Please ensure'...
             ' ttype is set correctly. Options are "WRAP" for a wrapping '...
             'transform and "ME" for a mirror-extended transform']);
          assert( isscalar(nbscales) && isscalar(nbangles),['Please ensure'...
             ' nbscales and nbangles are scalar values']);
		  % assert( subscale <= nbscales || subscale > 0 ,['Please ensure'...
		  % 	             ' subscale is a positive number and less that nbscales']);

          % Compute length of curvelet coefficient vector
          if strcmp(ttype,'ME')
             C = mefcv2(randn(m,n),m,n,nbscales,nbangles);

             hdr{1}{1} = size(C{1}{1});
             cn = prod(hdr{1}{1});
             for i = 2:nbscales
                nw = length(C{i});
                hdr{i}{1} = size(C{i}{1});
                hdr{i}{2} = size(C{i}{nw/2+1});
                cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
             end
          else
             [tmphdr, cn] = fdct_sizes_mex(m,n,nbscales,nbangles,logical(finest));
             hdr = cell(1,nbscales);
             hdr{1} = {[tmphdr{1:2}]}; 
             for i = 2:nbscales - (~finest)
                j = 3 + 5*(i-2);
                hdr{i}={[tmphdr{j+1:j+2}];[tmphdr{j+3:j+4}];[tmphdr{j}]};
             end
             if ~finest,  hdr{end} = {[tmphdr{end-1:end}];1};       end;
          end

          % Construct operator
          op = op@opSpot('Scale', cn, cn);
          op.cflag     = 1;
          op.nbscales = nbscales;
          op.nbangles = nbangles;
          op.finest = finest;
		  op.subscale = subscale;
          op.header = hdr;
          op.nbcoeffs = cn;
          op.dims = [m,n];
          op.ttype = ttype;
       end % Constructor

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % rrandn             
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % overloaded to produce a vector that really falls in the range of op
       function y = rrandn(op)
          y = op.drandn;
          y = multiply(op,y,1);
       end
       
    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function x = multiply(op,x,mode)
	   x_temp = zeros(size(x));
         if mode == 1 
            if strcmp(op.ttype,'ME')
               x = spot.utils.mefdct_v2c(x,op.header,op.nbangles);
               x_temp = spot.utils.mefdct_v2c(x_temp,op.header,op.nbangles);
            else
               x = spot.utils.fdct_v2c(x,op.header);
               x_temp = spot.utils.fdct_v2c(x_temp,op.header);
            end

            for mm = op.subscale
			   for nn = 1:length(x{mm})
				   x_temp{mm}{nn} = x{mm}{nn};
			   end
			end
	        x = spot.utils.fdct_c2v(x_temp,op.nbcoeffs);
         else 
            if strcmp(op.ttype,'ME')
               x = spot.utils.mefdct_v2c(x,op.header,op.nbangles);
               x_temp = spot.utils.mefdct_v2c(x_temp,op.header,op.nbangles);
            else
               x = spot.utils.fdct_v2c(x,op.header);
               x_temp = spot.utils.fdct_v2c(x_temp,op.header);
            end

            for mm = op.subscale
			   for nn = 1:length(x{mm})
				   x_temp{mm}{nn} = x{mm}{nn};
			   end
			end
	        x = spot.utils.fdct_c2v(x_temp,op.nbcoeffs);
         end
         
       end % Multiply
       
    end % Methods
   
end % Classdef