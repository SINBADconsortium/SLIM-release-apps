function op = opDataMatrix_MS(A)

nr = size(A,2);
ns = size(A,3);
nf = size(A,1);

AP = permute(A, [2 3 1]);
clear A

subfunc_handle = @(x,mode) opDataMatrix_intrnl(x,mode);
m = nr*ns*nf;
n = nr*nr*nf;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction

	% Subfunction
	function y = opDataMatrix_intrnl(x,mode)

		switch mode  
		    case 0
		        c =~isreal(AP);
		        y = {m,n,[0,1,0,1],opinfo};
    
		    case 1
		        x = reshape(x,nf,nr,nr);
		        x = permute(x,[2 3 1]);
  
		        y = zeros(nr,ns,nf);
  
		        for k=1:nf
		            y(:,:,k) = (x(:,:,k)*AP(:,:,k));
		        end
  
		        y = permute(y,[3 1 2]);
		        y = y(:);
    
		    case 2
		        x = reshape(x,nf,nr,ns);
		        x = permute(x,[2 3 1]);

		        y = zeros(nr,nr,nf);

		        for k=1:nf
		            y(:,:,k) = (x(:,:,k)*(AP(:,:,k)'));
		        end

		        y = permute(y,[3 1 2]);
		        y = y(:);
		end
	end
end