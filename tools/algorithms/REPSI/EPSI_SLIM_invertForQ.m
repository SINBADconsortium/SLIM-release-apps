function [update, update_unwindowed] = EPSI_SLIM_invertForQ(b, opEPSI, x_phys, nt, nr, ns, nt_conv, dt, q, wavelet_window_start, wavelet_window_end, wavelet_reg_lambda, options)
% Match for wavelet on residual data b
    
    % re-shape vec-ed datacube to physical dimensions
    if options.parallel
        x_phys = distVec2distArray(x_phys, [nr, ns, nt]);
    else
        x_phys = reshape(x_phys, [nt, nr, ns]);
    end
    
    % define the operator for least-squares matching
    [temp1 temp2 temp3 OP_CONV_WITH_DATA] = opEPSI(x_phys, q, 0, dt, 0);


    % garbage collection
    clear temp1
    clear temp2
    clear temp3
    clear x_phys


    % calculate mastering window to eliminate everything outside a predifined source wavelet time-length
    taper_length = floor(options.wavelet_taperperc * floor(options.q_estLength_posT ./ dt));
    if (taper_length <= 0)
        % use hard thresholding window
        OP_S_WIND = opTimeWindow(zeros(size(OP_CONV_WITH_DATA,2),1), wavelet_window_start, wavelet_window_end);
    else
        % use a window with cosine tapers       
        OP_S_WIND = opTimeWindow_Tapered(zeros(size(OP_CONV_WITH_DATA,2),1), wavelet_window_start, wavelet_window_end, taper_length);
    end
    
    % The system to invert is OP_CONV_WITH_DATA
    A = OP_CONV_WITH_DATA;
    
    % wavelet matching in the time domain via LSQR, obtain approximate pseudo-inverse of Ax = b using LSQR
    % where A defines the operation of convolving a source signature x with every trace of x_phys (current estimate for G)
    % and b is the data-space residual of previous estimates
    A_adaptor = func_adaptor_lsqr(A); % (LSQR is incompatible with native SPOT operators, use an adaptor)
    update_unwindowed = lsqrms(size(A,1),size(A,2),A_adaptor,0,0,b,wavelet_reg_lambda,0,0,0,5,options.verbosity);
    
    % zero-out everything outside of the predifined source wavelet time-length
    update = OP_S_WIND * update_unwindowed(:);
    
    % exact line search scaling of the solution
    x = update;
    
    g = A' * b;
    y = A * x;    
    s = (x' * g) / (y' * y);
    
    s = abs(undist(s)); disp(['   ... sacled by: ' num2str(s)])
    update = x * s;

end