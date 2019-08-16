function [OP_EPSI OP_PTerm OP_QTerm] = opEPSI_MS_dist(data,q,dt,use_oblique,src_mask,F_taper)
% opEPSI for missing source acquisition
% no topmute or time-padding, consistent with MIGSRM
% alpha: alpha for Tukeywin in MIGSRM, default 0 here
    
    if not(exist('use_oblique','var'))
        use_oblique = 0;
    end
    
    %% Get relavant dimension information (assumes data is: d3=Time, d1=reciever, d2=shots)
    dims = size(data);
    nt   = dims(3);
    nr   = dims(1);
    ns   = dims(2);

    %% Make specialized FFT operators for the convolution and data transformations
    
    % determine number of frequencies 
    F1D = opFFTsym_conv(nt);
    nf = size(F1D,1);
    
    % construct the DFT operator that acts on the whole datacube for convolution
    F_P = oppKron2Lo(F1D,opDirac(nr*ns));
    F_G = oppKron2Lo(F1D,opDirac(nr*nr));

    %% apply correction for obliquity factor in 2D data
    if use_oblique
        funcTimeConv = convolution_time2D(nr,ns);
        obliquity_factors = make_obliquityFactor(dt,nf);
        obliquity_factors = distributed(obliquity_factors.');
        OPOBLIQ = oppDistFun(obliquity_factors,funcTimeConv,0);
    	OPOBLIQINV = opObliqInv(dt,nf,1);
    else
        OPOBLIQ = opDirac(nf*nr*ns);
        OPOBLIQINV = opDirac(nf);
    end
    
    %% Prepare data used for the multidimentional data-data convolutioon
    
    % Define the mastering window the the finite Fourier series convolution (using a Tukey consine window)
    % Time domain windowing
    time_taper_frac = 0.1; % ranges as a fraction from 0 to 1
    window_start = 1;
    window_end = nt;
    window_length = window_end - window_start + 1;
    taper_length = floor((time_taper_frac/2) * window_length);
    OPWIND_TIME = opTimeWindow_Tapered(zeros(nt,1),window_start,window_end,taper_length);

    % Receiver coordiante windowing
    recv_taper_frac = 0.05; % ranges as a fraction from 0 to 1
    window_start = 1;
    window_end = nr;
    window_length = window_end - window_start + 1;
    taper_length = floor((recv_taper_frac/2) * window_length);
    OPWIND_RECV = opTimeWindow_Tapered(zeros(nr,1),window_start,window_end,taper_length);
    
    % Window the data
    OPWIND = oppKron2Lo(OPWIND_TIME, opKron(opDirac(ns),OPWIND_RECV));
    data = OPWIND*data(:);
    
    % Make frequency slices for convolution
    P_f = oppKron2Lo(opDiag(F_taper), opDirac(nr*ns))*F_P*data(:);
    clear data
    P_f = OPOBLIQ*P_f;
    P_f = reshape(P_f,nr,ns,nf);
    
    %% Define data for the wavelet term Q
    q_f = opDiag(F_taper)*F1D*q;
    
    %% P and Q operator
    funcMultidimConv = MDconv_2D_MS(nr,ns);
    P = oppDistFun(P_f, funcMultidimConv, 0);
    Q = oppKron2Lo(opDiag(q_f),opKron(opRestriction(nr,src_mask),opDirac(nr)));
    
    %% Combine (Q-P) and a muting operator to make the EPSI upgoing wavefield prediction operator
    OP_EPSI = F_P'*(Q-P)*F_G;
    OP_PTerm = -F_P'*P*F_G;
    OP_QTerm = F_P'*Q*F_G;