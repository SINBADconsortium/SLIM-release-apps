function op = opEnsureReal(n)

% OPENSUREREAL    Returns an operator that acts like an Identity operator, but
%                 ensures that the output is always real. Used to clean-up small
%                 imaginary parts of physical signals that result from
%                 computational rounding errors in spectrum analysis/operations
	
subfunc_handle = @(x,mode) opEnsureReal_intrnl(x,mode);

    op = opFunction(n,n,subfunc_handle); % return a SPOT operator using constructor opFunction

    function y = opEnsureReal_intrnl(x,mode)
    	if mode == 0
    	        y = {n,n,[1,1,1,1],{'opEnsureReal'}};
    	else
    	        y = real(x);
    	end
    end
end

