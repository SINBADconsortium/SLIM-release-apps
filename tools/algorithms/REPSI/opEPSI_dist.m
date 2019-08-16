function [OP_EPSI, OP_P_TERM, OP_Q_TERM, OP_CONV_WITH_DATA, OP_OBLIQ, OP_OBLIQINV_Q] = opEPSI_dist(data,q,topmute,dt,use_oblique)
%% OPEPSI      The EPSI observation operator on 2D survey datacube (takes the data always in time domain)
%  usage: opEPSI(data,s,topmute,dt,use_oblique)
%  
% DISTRIBUTED VERSION: requires "data" to be in a distributed 3D cube with time on the last axis (time-slice gathers),
%                      which should also be the axis that the data is distributed over.
% 
% PARAMETERS:
%     USE_OBLIQUE         1 to include obliquity factor int he data (phase shift), 0 to disable


	if ~exist('use_oblique','var')
		use_oblique = 0;
	end
	
	q = undist(q);
	
	%% Get relavant dimension information (assumes data is: d1=reciever, d2=shots, d3=time)
	dims = size(data);
	nr   = dims(1);
	ns   = dims(2);
	nt   = dims(3);
	nt_conv = 2*nt; % time sample length of padded kernel for convolution
    
    data = distVectorize(data);
    
	%% Make specialized FFT operators for the convolution and data transformations
	
	% determine number of frequencies 
    F1D = opFFTsym_conv(nt_conv);
	nf = size(F1D,1);

	%% Padding and chopping operators for non-wrap-around Fourier domain convolution
	OPPAD_BOTTOM = opPadBottom([nt 1],nt_conv,1);
	OPPAD_TOP = opPadTop([nt 1],nt_conv,1);
	OPCHOP_TOP = OPPAD_TOP';
    OPCHOP_BOTTOM = OPPAD_BOTTOM';
    
    % construct the DFT operators that acts on the whole datacube for (non-circular) multidimensional convolution
	F = oppKron2Lo( F1D * OPPAD_BOTTOM , opDirac(nr*ns) );
	invF = oppKron2Lo( OPCHOP_TOP * F1D' , opDirac(nr*ns) ); % for opFFTsym_conv the adjoint mode implements the inverse
    
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

	%% Define an operator to block event near t=0, to prevent inverting energy to acquisition source (trivial solution)
    OPMUTE = oppKron2Lo(opMask(nt,[topmute+1:nt]),opDirac(nr*ns));

	%% Prepare data used for the multidimentional data-data convolution
	
    % Define the mastering window the the finite Fourier series convolution (using a Tukey consine window)
        % Time domain windowing
        time_taper_frac = 0.1; % ranges as a fraction from 0 to 1
    	window_start = ceil((topmute+1) * 0.8);
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
    data = OPWIND * data;
        
    % Make frequency slices of the data for convolution
    F_data = oppKron2Lo(F1D * OPPAD_TOP, opDirac(nr*ns));
	data_conv_f = F_data * data; % Convolution is in Fourier domain
	clear data
	
	% An obliquity factor may be needed so that P = X0*SI - OBLIQ*X0*P.
	data_conv_f = OPOBLIQ * data_conv_f;
	data_conv_f = distVec2distArray(data_conv_f,[nr ns nf]);	 
    
	%% Define data for the wavelet term Q
	if length(q) < nt % make sure s is the right length
		q = [q(:); zeros(nt-length(q),1)];
	end
	OPPAD_WAVELET = opPadTop(q,nt_conv);
	
	if length(q) == nt_conv
	    qpad = q; % q alrady has both causal and anti-causal parts
    else
    	qpad = OPPAD_WAVELET * q; % wavelet causal, pad zeros to anti-causal part
	end
	
	qf = F1D * qpad; % Fourier coefficients
	qf = distributed(qf.');
	
	%% Finally the multidimentional data-data convolution operators can be constructed
	funcMultidimConv = convolution_multidim2D(nr,ns);
	funcTimeConv = convolution_time2D(nr,ns);
	funcTimeConvWithData = convolutionWithData_time2D(nr,ns);
	
	P = oppDistFun(data_conv_f, funcMultidimConv, 0);
	Q = oppDistFun(qf, funcTimeConv, 0);
	
	
	%% Combine (Q-P) and a muting operator to make the EPSI upgoing wavefield prediction operator
    OP_EPSI = invF * (Q - P) * F * OPMUTE;

    
    if nargout > 1
        % additional operators used for matching
        OP_P_TERM = invF * P * F;
        OP_Q_TERM = invF * Q * F;
        OP_CONV_WITH_DATA = oppKron2Lo( OPCHOP_BOTTOM * F1D' , opDirac(nr*ns) ) * oppDistFun(data_conv_f, funcTimeConvWithData, 3) * F1D;
        OP_OBLIQ = F' * OPOBLIQ * F;
        OP_OBLIQINV_Q = F1D' * OPOBLIQINV * F1D;
    end

end