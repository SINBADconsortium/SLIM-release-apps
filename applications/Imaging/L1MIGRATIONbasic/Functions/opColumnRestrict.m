classdef opColumnRestrict < opSpot
% OPCOLUMNRESTRICT   Restriction operator on matrix columns
%
%    OPCOLUMNRESTRICT(M,N,IDX,TYPE), with TYPE = 'DISCARD', creates
%    an operator that extracts the columns indicated by IDX from a 
%    given M by N matrix. The adjoint operator takes an M by
%    length(IDX) matrix and outputs an M by N matrix with the
%    columns filled by the given input matrix. Note that all input
%    and output matrices are in vectorized from. When TYPE =
%    'ZERO', all columns that are not in IDX are zero-padded
%    instead of discarded.
%
%    See also opMask, opRestriction.


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
       function op = opColumnRestrict(m,n,idx,type)

		if (min(idx) < 1) || (max(idx) > n)
		   error('Index parameter must be integer and match dimensions of the data');
		end

		if (nargin < 4) || isempty(type)
		  type = 'discard';
		end

		if strcmp(lower(type),'discard')
		   l  = length(idx);
		   fun = @(x,mode) opColumnRestrict_intrnl1(idx,l,m,n,x,mode);
		  % Construct operator
          idex1 = m*l; idex2 = m*n;
		else
		   invidx = ones(n,1);
		   invidx(idx) = 0;
		   idx = find(invidx);
		   fun  = @(x,mode) opColumnRestrict_intrnl2(idx,m,n,x,mode);
           idex1 = m*n; idex2 = m*n;
		end
		
		  % Construct operator
		
	      op = op@opSpot('ColumnRestrict', idex1, idex2);
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


function y = opColumnRestrict_intrnl1(idx,l,m,n,x,mode)

if mode == 1
   y = reshape(x,m,n);
   y = reshape(y(:,idx),m*l,1);
else
   y = zeros(m,n);
   y(:,idx) = reshape(x,m,l);
   y = reshape(y,m*n,1);
end
end
function y = opColumnRestrict_intrnl2(invidx,m,n,x,mode)

   y = reshape(x,m,n);
   y(:,invidx) = 0;
   y = y(:);
end