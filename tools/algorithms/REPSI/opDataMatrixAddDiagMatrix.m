function op = opDataMatrixAddDiagMatrix(A,diag_coefs,scale,opinfo)
% A fast implementation of the following SPARCO lines using functions here:
%
%   A = opDataMatrix(A,scale);
%   B = opDiag3D(size(A,1),size(A,2),size(A,3),diag_coefs);
%   OP = opSum(A,B, [scale,1])
%
%   scale defaults to 1

%  assertions: diag_coefs is a vector of length size(A,1) or a scalar
%  A is a 3D array

% metainfo, for record-keeping
if exist('opinfo') ~= 1
  opinfo = {'DataMatrixAddDiagMatrix', []};
elseif isstr(opinfo)
  opinfo = {'DataMatrixAddDiagMatrix', opinfo};
end



% defaults on scale
if exist('scale') ~= 1
    scale = 1;
end

nr  = size(A,2);
ns  = size(A,3);
nf = size(A,1);

% defaults for diag_coefs
if isscalar(diag_coefs)
    diag_coefs = diag_coefs * ones(nf,1);
else
    assert(length(diag_coefs) == nf);
    diag_coefs = diag_coefs(:);
end

% data-reshape, and force execution by clearing A
AP = permute(A, [2 3 1]);
clear A;

subfunc_handle = @(x,mode) opDataMatrix_intrnl(x,mode);
m = nr*ns*nf;
n = nr*ns*nf;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction
	% Subfunction
	function y = opDataMatrix_intrnl(x,mode)

		switch mode  
		    case 0
		        c =~isreal(AP);
		        y = {m,n,[0,1,0,1],opinfo};
    
		    case 1
		        x = reshape(x,nf,nr,ns);
		        x = permute(x,[2 3 1]);
  
		        for k=1:nf
		            x(:,:,k) = (diag_coefs(k) .* x(:,:,k)) + (scale .* (x(:,:,k)*AP(:,:,k)));
		        end
  
		        x = permute(x,[3 1 2]);
		        x = x(:);
		        y = x;
    
		    case 2
		        x = reshape(x,nf,nr,ns);
		        x = permute(x,[2 3 1]);

		        for k=1:nf
		            x(:,:,k) = (conj(diag_coefs(k)) .* x(:,:,k)) + (conj(scale) .* (x(:,:,k)*(AP(:,:,k)')));
		        end

		        x = permute(x,[3 1 2]);
		        x = x(:);
		        y = x;
		end
	end
end