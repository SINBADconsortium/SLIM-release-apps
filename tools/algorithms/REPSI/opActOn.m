function op = opActOn(n1,n2,n3,data,scaling)
% OPACTON         Hadamard-product of each element of x to the corresponding 2D-slice along 1st dimension of a 3D data
%                 (for the source signature term in EPSI_SLIM)

data = reshape(data,n1,n2,n3);
data = permute(data, [2 3 1]);

if exist('scaling') ~= 1
    scaling = 1;
end


subfunc_handle = @(x,mode) opActOn_intrnl(x,mode);
m = n1*n2*n3;
n = n1;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction

	function y = opActOn_intrnl(x,mode)

		switch mode
		    case 0
		        y = {m,n,[0,1,0,1],{'opActOn'}};
   
		    case 1
		        x = x .* scaling;
		        x = repmat(reshape(x, [1 1 n1]), [n2 n3]);
		        x = x .* data;
		        x = permute(x, [3 1 2]);
		        y = x(:);

		    case 2
		        x = x .* conj(scaling);
		        x = reshape(x,n1,n2,n3);
		        x = permute(x, [2 3 1]); 
		        x = x .* conj(data);
		        x = reshape(x,[],n1);
		        y = sum(x,1);
		        y = y(:);        
		end
	end
end
