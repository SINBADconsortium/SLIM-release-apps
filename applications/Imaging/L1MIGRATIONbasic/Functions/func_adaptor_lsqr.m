function op = func_adaptor_lsqr( fname )
% Adapts a SPARCO operator calling to the format expected by LSQR

op = @(x,mode) spotadaptor(x,mode,fname);


end


function y = spotadaptor(x,mode,fname);
if mode == 1
    y = fname*x;
elseif mode == 2
    y = fname'*x;
end
end

