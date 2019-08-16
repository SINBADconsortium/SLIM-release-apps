classdef opWaveAtom < opSpot
%OPWAVEATOM  Two-dimensional curvelet operator.
%
%   Spot operator warp for 2D forward mirror-extended wave atom transform
%   Forward mode is Forward wave atom transfrom; while adjoint mode is 
%	inverse transform
%
%   Input arguements:
%   	N:   size of the images which should be square, while N should be
%   	     a power of 2
%   	ttyple: whether use mirror-extended
%       	'NO': Normal wave atom transform, not mirror-extended
%	    	'ME': mirror-extended wave atom transform
%   	pat: specifies the type of frequency partition which satsifies
%   		 parabolic scaling relationship. pat can either be 'p' or 'q'.
%   	tp:  type of tranform.
%   		'ortho': frame based on the orthobasis construction of 
%   			the standard wave atom
%   		'directional': real-valued frame with single oscillation direction
%   		'complex': complex-valued frame
%
%   See also Waveatom.
%   auther: Xiang Li

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
	   nbcoeffs;
	   header;
	   ttype;
	   dims;
       pat;
       tp;
  
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opWaveAtom(n,ttype,pat,tp)

          assert( isscalar(n),['Please ensure'...
            ' sizes are scalar values']);
          if nargin < 2, ttype = 'NO';                               end;
          if nargin < 3, pat   = 'p';                                end;
          if nargin < 4, tp    = 'directional';                      end;
          

          assert( strcmp(ttype,'NO') || strcmp(ttype,'ME'), ['Please ensure'...
             ' ttype is set correctly. Options are "NO" for a normal'...
             'transform and "ME" for a mirror-extended transform']);

          % Compute length of curvelet coefficient vector
          if strcmp(ttype,'ME')
             C = mefwa2(zeros(n),pat,tp);
		  else 
			 C = fwa2(zeros(n),pat,tp);
          end

		  cn   = 0;
		  for i = 1:prod(size(C));
	      	 for j = 1:prod(size(C{i}))
				hdr{i}{j} = size(C{i}{j});
				
				cn  = cn + prod(size(C{i}{j}));
		   	 end
			 hdr{i} = reshape(hdr{i},size(C{i}));
		  end
		  
		  hdr = reshape(hdr,size(C));
          % Construct operator
          op = op@opSpot('WaveAtom', cn, n*n);

          op.header    = hdr;
          op.nbcoeffs  = cn;
          op.dims      = [n,n];
		  op.ttype     = ttype;
		  op.pat       = pat;
		  op.tp        = tp;

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
         if mode == 1
            % Analysis mode
			if strcmp(op.ttype,'ME')
				x = mefwa2(reshape(x,op.dims),op.pat,op.tp);
			else 
				x = fwa2(reshape(x,op.dims),op.pat,op.tp);
			end
			x = wa_c2v(x,op.nbcoeffs);
			
         else
	        % Synthesis mode 
			x = wa_v2c(x,op.header);
			if strcmp(op.ttype,'ME')
				x = meiwa2(x,op.pat,op.tp);
			else 
				x = iwa2(x,op.pat,op.tp);
			end
            x = x(:);
         end
         
       end % Multiply
       
    end % Methods
   
end % Classdef



%=== subfunction make a cell into a vector ==== 
function y = wa_c2v(x,cn)
	y = [];
	for m = 1:prod(size(x))
		for n = 1:prod(size(x{m}))
			s = x{m}{n};
			y = [y;s(:)];
		end
	end
	if length(y) ~= cn, warning('length of coeffiences does not match');end
end

%=== subfunction make a vec to be a cell aboarding to hdr==== 
function y = wa_v2c(x,hdr);
	for m = 1:prod(size(hdr));
		for n = 1:prod(size(hdr{m}));
			y{m}{n} = reshape(x(1:prod(hdr{m}{n})),hdr{m}{n}(1),hdr{m}{n}(2));
			x(1:prod(hdr{m}{n})) = [];
		end
		y{m} = reshape(y{m},size(hdr{m}));
	end
	y = reshape(y,size(hdr));
end