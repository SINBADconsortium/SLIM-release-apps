function subfunc_handle = convolution_time2D(nr, ns)
% Y = multidimConvolution2D(nr, ns)
% 
% (Acts on monochromatic slices of seismic line, or 2D, data)
% Returns a function for oppDistFun that implements the multi-dimensional
% convolution used in EPSI modeling operator. For split-spread data in A 
% (monochromatic according to detail-hiding operator) and primary IR in
% X, this implements Y = qX when mode = 1 and its adjoint operation when
% in mode = 2 (on X)
% 
% Y, X are matrices of size nr-by-ns, q is a scalar

    
subfunc_handle = @(q, x, mode) convolution_time2D_intrnl(q, x, mode);

    function y = convolution_time2D_intrnl(q, x, mode);

    switch mode
        case 0
            % [m n cflag linflag]
            y = [nr*ns, nr*ns, 1, 1];

        case 1
            x = q .* x;
            y = x;

        case 2
            x = conj(q) .* x;
            y = x;
    end
    end
end