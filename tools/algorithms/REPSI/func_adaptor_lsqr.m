function op = func_adaptor_lsqr( OP )
% Adapts a SPOT operator calling to the format expected by LSQR

op = @(mode,m,n,x,iw,rw) helper_func(x,mode);

function y = helper_func(x,mode)
    if mode == 1
        y = OP * x;
    elseif mode == 2
        y = OP' * x;
    end
end
end