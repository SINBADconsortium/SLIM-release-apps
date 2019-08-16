function op = opDataMatrix(A,scale,opinfo)
% OPDATAMATRIX      Used for Sparse Primary Inversion to convolve primary with data
%                   Right matrix-multiplies the data with the primary, as stacks of 2D slices
%                   on a per-slice basis

if exist('opinfo') ~= 1
  opinfo = {'DataMatrix', []};
elseif isstr(opinfo)
  opinfo = {'DataMatrix', opinfo};
end

if exist('scale') ~= 1
    scale = 1;
end

nr  = size(A,2);
ns  = size(A,3);
nf = size(A,1);

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
  
		        y = zeros(nr,ns,nf);
  
		        for k=1:nf
		            y(:,:,k) = scale.* (x(:,:,k)*AP(:,:,k));
		        end
  
		        y = permute(y,[3 1 2]);
		        y = y(:);
    
		    case 2
		        x = reshape(x,nf,nr,ns);
		        x = permute(x,[2 3 1]);

		        y = zeros(nr,ns,nf);

		        for k=1:nf
		            y(:,:,k) = conj(scale).* (x(:,:,k)*(AP(:,:,k)'));
		        end

		        y = permute(y,[3 1 2]);
		        y = y(:);
		end
	end
end