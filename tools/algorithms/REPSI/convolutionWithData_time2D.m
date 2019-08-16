function subfunc_handle = convolutionWithData_time2D(nr, ns)
% Y = convolutionWithData_time2D(nr, ns)
% 
% (Acts on monochromatic slices of seismic line, or 2D, data)
% Returns a function for oppDistFun that implements a time-convolution 
% with data used in EPSI modeling operator. For split-spread data in A 
% (monochromatic according to detail-hiding operator) and source wavelet in
% x, this implements Y = A.*x when mode = 1 and its adjoint operation when
% in mode = 2 (on x)
% 
% Y, A are matrices of size nr-by-ns, x is a scalar

    
subfunc_handle = @(A, x, mode) convolutionWithData_time2D_intrnl(A, x, mode);

    function y = convolutionWithData_time2D_intrnl(A, x, mode);

    switch mode
        case 0
            % [m n cflag linflag]
            y = [nr*ns, 1, 1, 1];

        case 1
            y = x .* A;
            y = y(:);

        case 2
            x = x .* conj(A(:));
            y = sum(x);
    end
    end
end