function [OP_EPSI, OP_P_TERM, OP_Q_TERM, OP_CONV_WITH_DATA, OP_OBLIQ, OP_OBLIQINV_Q] = opEPSI_serial(data,q,topmute,dt,use_oblique)
%% OPEPSI      The EPSI observation operator on 2D survey datacube (takes the data always in time domain)
%  usage: opEPSI(data,s,topmute,dt,use_oblique)
%  
% PARAMETERS:
%     USE_OBLIQUE         1 to include obliquity factor int he data (phase shift), 0 to disable


	if ~exist('use_oblique','var')
		use_oblique = 0;
	end
	
	%% Get relavant dimension information (assumes data is: d1=Time, d2=reciever, d3=shots)
	dims = size(data);
	nt   = dims(1);
	nr   = dims(2);
	ns   = dims(3);
	nt_conv = 2*nt; % time sample length of padded kernel for convolution

	%% Make specialized FFT operators for the convolution and data transformations
	
	% determine number of frequencies 
    F1D = opFFTsym_conv(nt_conv);
	nf = size(F1D,1);
	
	% construct the DFT operator that acts on the whole datacube for convolution
	F = opFFTsym_conv_datacube(nt_conv,nr,ns);
	invF = F'; % for opFFTsym_conv_datacube the adjoint mode implements the inverse

	%% Padding and chopping operators for non-wrap-around Fourier domain convolution
	OPPAD_BOTTOM = opPadBottom(data,nt_conv);
	OPPAD_TOP = opPadTop(data,nt_conv);
	OPCHOP_TOP = OPPAD_TOP';
    OPCHOP_BOTTOM = OPPAD_BOTTOM';
    
    %% apply correction for obliquity factor in 2D data
	if use_oblique
    	OPOBLIQ = opObliq(dt,nf,nr*ns);
    	OPOBLIQINV = opObliqInv(dt,nf,1);
	else
	    OPOBLIQ = opDirac(nf*nr*ns);
	    OPOBLIQINV = opDirac(nf);
	end
    
	%% Define an operator to block event near t=0, to prevent inverting energy to acquisition source (trivial solution)
    OPMUTE = opTriMute(data,topmute,0,0);

	%% Prepare data used for the multidimentional data-data convolutioon
    
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
    OPWIND = opKron(opDirac(ns),OPWIND_RECV,OPWIND_TIME);
    data = OPWIND * data(:);
    data = reshape(data,nt,nr,ns);
    
    % Make padded versions of the data for convolution
	data_conv = OPPAD_TOP * data(:);
	clear data
	data_conv_f = F * data_conv(:); % Convolution is in Fourier domain
	clear data_conv
	% An obliquity factor may be needed so that P = X0*SI - OBLIQ*X0*P.
	data_conv_f = OPOBLIQ * data_conv_f;
	data_conv_f = reshape(data_conv_f,nf,nr,ns);	 

    
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
	
	
	%% Finally the multidimentional data-data convolution operators can be constructed
    QminusP = opDataMatrixAddDiagMatrix(data_conv_f,qf,-1);
    
        % These are only needed if intermediate operators are needed
        if nargout > 1
            P = opDataMatrix(data_conv_f);
        	Q = opDiag3D(nf,nr,ns,qf);
        end
    
	%% Combine (Q-P) and a muting operator to make the EPSI upgoing wavefield prediction operator
    OP_EPSI = OPCHOP_TOP * invF * QminusP * F * OPPAD_BOTTOM * OPMUTE;
    
    if nargout > 1
        % additional operators used for matching
        OP_P_TERM = OPCHOP_TOP * invF * P * F * OPPAD_BOTTOM;
        OP_Q_TERM = OPCHOP_TOP * invF * Q * F * OPPAD_BOTTOM;
        OP_CONV_WITH_DATA = OPCHOP_BOTTOM * invF * opActOn(nf,nr,ns,data_conv_f(:),1) * F1D;
        OP_OBLIQ = OPCHOP_BOTTOM * invF * OPOBLIQ * F * OPPAD_BOTTOM;
        OP_OBLIQINV_Q = F1D' * OPOBLIQINV * F1D;
    end

	clear data_conv_f

end