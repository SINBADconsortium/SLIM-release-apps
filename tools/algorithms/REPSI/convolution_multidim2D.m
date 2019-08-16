function subfunc_handle = convolution_multidim2D(nr, ns)
% Y = convolution_multidim2D(nr, ns)
% 
% (Acts on monochromatic slices of seismic line, or 2D, data)
% Returns a function for oppDistFun that implements the multi-dimensional
% convolution used in EPSI modeling operator. For split-spread data in A 
% (monochromatic according to detail-hiding operator) and primary IR in
% X, this implements Y = X*A when mode = 1 and its adjoint operation when
% in mode = 2 (on X)
% 
% Y, A, X are matrices of size nr-by-ns (for split spread require nr = ns)

    
subfunc_handle = @(A, x, mode) convolution_multidim2D_intrnl(A, x, mode);

    function y = convolution_multidim2D_intrnl(A, x, mode);

    switch mode
        case 0
            % [m n cflag linflag]
            y = [nr*ns, nr*ns, 1, 1];

        case 1
            x = reshape(x, [nr ns]);
            x = x * A;
            x = x(:);
            y = x;

        case 2
            x = reshape(x, [nr ns]);
            x = x * A';
            x = x(:);
            y = x;
    end
    end
end