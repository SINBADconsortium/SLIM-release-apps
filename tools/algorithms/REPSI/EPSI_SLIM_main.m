function EPSI_SLIM_main(varargin)

% Robust primary estimation by sparse inversion (rEPSI) using L1 minimization as sparsity regularization (via SPGL1)

% Author      : Tim Lin
%               Seismic Laboratory for Imaging and Modelling
%               Department of Earth & Ocean Sciences
%               The University of British Columbia
%         
% Date        : Faburary, 2011
% 
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% Requirements: SPGL1 (http://www.cs.ubc.ca/labs/scl/index.php/Main/Spgl1)
%               RSF (http://rsf.sourceforge.net/) with Matlab API (or ability for RSF file IO)
%               Sparco (http://www.cs.ubc.ca/labs/scl/sparco/)


% Test to see if using MATLAB R2009b or later
if verLessThan('matlab','7.9.0')
    error('EPSI_SLIM:MatlabVersionError', 'Your version of MATLAB is too old. EPSI_SLIM requires MATLAB version 7.9 (2009b) or later.')
end

% Parse options and specify their defaults, best read with syntax-highlighitng editor and line-wrapping turned OFF
[options] = process_options_to_struct(varargin, ...
...
...     % Option Name           % Default value
...
... % REQUIRED OPTIONS ========================
...     % Input data (see README regarding requirements on input data)
        'input_file',           '../data/Synthetic_long.mat', ...    % name of seismic datafile to process, can be .mat, .su (native), or .bin (presume axis1=time, axis2=reciever, axis3=shot)
...         % (the following options are mainly only needed for troubleshooting the input file)
        'input_ns',             [], ...         % number of time samples per trace, this is required to read from a pure binary datafile
        'input_endian',         [], ...         % can specify the endianess of binary/su datasets if different from native, 'b' for big, 'l' for little
        'dt',                   [], ...         % expected time-sample rate of data in seconds (overrides dt read from the data file)
        'test_readData',        0, ...          % if set to 1, launch into tesing mode to read data and display example gathers, then quit
            ...
...     % Required data-dependent information
        'topmuteT',             0.10, ...       % global muting time before first reflected arrival, setting this close to the first arrival prevents trivial solution
        'q_estLength_posT',     0.08, ...       % expected length of source signature wavelet in seconds (in positive time direction)
        'q_estLength_negT',     0.08, ...       % expected length of source signature wavelet in seconds (in negative time direction)
            ...
...     % Operation mode
        'useSparsity',          0, ...          % set to 1 if the solution should be seeked in a sparsifying representation (hybrid Wavelet and Curvelet domain), uses a lot of memory unless 'parallel' is used
        'parallel',             0, ...          % set to 1 if MATLAB parallel toolbox workers are available and parallel pool is active. Will execute in multi-process mode.
            ...
...     % Output files, if empty then no files will be output
        'output_primary_file',  [], ...         % result estimated primary wavefield, can be .mat, .su (native), or .bin
        'output_primaryIR_file',[], ...         % result estimated primary impulse response (uncolvolved with src sig), can be .mat, .su (native), or .bin
        'output_wavelet_file',  [], ...         % result estimated global source signature wavelet, can be .mat, .su (native), or .bin
...
... % OPTIONAL PARAMS ========================
...     % Tweaks and inversion parameters that alter result quality and algorithm operation
        'relError',             0.05, ...       % relative 2-norm tolerance level of residual (relative to 2-norm of the data), used for 'sigma'
        'maxTotalIter',         70, ...         % maximum number of total L1 gradient iterations to use (each gradient step costs 1 SRME-type multidimensional convolution)
        'maxInitialIter',       5, ...          % maximum number of initial L1 gradient iterations before initial spike-picking
            ...
        'initialSpikeAmp',      1.0, ...        % affects how aggressively anti-causal wavefields are ignored. Values >1 is better at multiple removal but may include non-causal artifacts. Normal range is 0.5~1.5
            ...
        'useOblique',           0, ...          % if =1, turn on obliquity factor (same role as in traditional SRME processing, to account for deghosted data)
            ...
        'invertQ',              1, ...          % turns on/off source signature estimation, if =0 then please make sure 'wavelet_file' is specified
        'wavelet_file',         [], ...         % optional estimation of a source signature file (currently only supports a mat file that stores this in a variable 'q')
            ...
        'wavelet_reg_lambda',   1, ...          % if >0, defines a damping parameter lambda on the 2-norm of the source signature when performing least-squares matching
        'wavelet_taperperc',    0.2, ...        % fraction of q_estLength_posT to taper-edge the estimates source by, if >0 then a cosine taper is used instead of hard threshold
        'window_startT',        0, ...          % if >0, indicate a time window where the initial L1 gradient iterations before initial spike-picking is done, in seconds
        'window_endT',          0, ...          %   (ideally 'window_startT' and 'window_endT', if specified, should roughtly window the first strongly reflected ocean bottom event)
            ...
        'polishWavelet',        1, ...          % if 1, match for the source signature one final time right before outputting results
            ...
        'minL1iters',           5, ...          % impose minimum number of L1 iterations within each alternating optimization loop
        'maxOuterIter',       inf, ...          % if <inf, puts limit on number of alternating optimization (block coordinate descent) loops
            ...
        'subspaceMin',          0, ...          % 1 enables the us eof subspace minimization step in SPGL1
            ...
        'curvType',        'WRAP', ...          % Curvelet type (when useSparsity=1), can be 'WRAP' (wrapping) or 'ME' (mirror-extended)
            ...
        'singlePrecision',      0, ...          % set to 1 to process everything in single precision (and output in single prec)       
            ...
            ...
...     % Optional preprocessing on the input data
        'downshift',            0, ...          % number of zero time sample padded to the top of the data
        'padtime',              0, ...          % number of zero time sample padded to the bottom of data
        'timeweightgamma',      0, ...          % if >0, time-weights data with exp( t * timeweightgamma ) during processing
        'renormData',           0, ...          % renormalize data such that maximum value is 1
            ...
        'num_traces_missing',   0, ...          % used to experiment with in-filling missing near offsets, not yet fully implemented
            ...
            ...
...     % Debug and developmental use
        'show_preview',         0, ...          % if =1, displays a graphical summary of results once processing is done
            ...
        'savepreviewmat',       0, ...          % save preview files in mat format
        'preview_file',         'testmodel_preview.mat', ... % target mat filename name of a single-shot preview of the solution
        'preview_diagonal',     0, ...          % 0 for 0 offset gathers, >0 for offset in dr
        'savesolmat',           0, ...          % if ~=0, saves all internal solution in mat format in sol_file
        'sol_file',             'testmodel_sol.mat', ...
        'savecoefmat',          0, ...          % if ~=0, saves solution coefficients in mat format in coef_file
        'coef_file',            'testmodel_coef.mat', ... % name of coef to save
        'savepickedspikes',     0, ...
        'verbosity',            1 ...           % affects internal solver verbosity
);


% add correct paths if ran interactively:
if nargin == 0
    addpath(genpath('.'))
end

% Check to see that SLIM version of SPGl1 is used, if not in data-reading mode
if (options.test_readData == 0)
    try
        assert(spgl1('is_SLIM_version') == 1, 'ERROR: Cannot find SLIM version of SPGl1, please re-examine the MATLAB look-up path to make sure it is included')
    catch ME
        error('ERROR: Cannot find SLIM version of SPGl1, please re-examine the MATLAB look-up path to make sure it is included')
    end
end

% load seismic data to process
[data file_dt] = read_EPSIdata(options.input_file, options.input_ns, options.input_endian);

% load source signature if defined (currently only works with .mat files)
if ~isempty(options.wavelet_file) % source defined
    load(options.wavelet_file,'q');
end

% determine dt
if ~isempty(options.dt) && (options.dt ~= 0)
    dt = options.dt;    
elseif ~isempty(file_dt) && (file_dt ~= 0)
    dt = file_dt;
else
    error('dt is never specified, neither in the data file nor in the options')
end

% check if this is a data-reading test
if (options.test_readData ~= 0)
    %display relavant info
    dt
    display_dataCube(data)
    return % quit StabalizedEPSI
end

%% Data pre-processing -----------------------
    
    % nt = n_timeSamples, ns = n_shots, nr = n_recievers
    [data, nt, ns, nr, options] = EPSI_SLIM_preprocessData(data, dt, options);
    
    % NOTE: if parallel processing is turned on, the data is also now sorted in (and distributed over) time-slice gathers

%% Initialize certain variables and data statistics -----------------------
    
    % nt_conv is the padded time domain used to carry out the non-circular data convolution
    nt_conv = 2*nt;
    
    % initialize s if undefined
    if isempty(options.wavelet_file)
        if options.invertQ
            q = zeros(nt,1);
        else
            error('Since "invertQ" is off, a known source signature is required');
        end
    end

    % make sure s is of the right length
    if length(q) < nt
        q = [q(:); zeros(nt-length(q),1)];
    end
    
    % total gradient iterations of the primaryIR matching step
    total_itercount = 0;
    
    data_1norm = norm(vectorizeData(data),1);
    data_infnorm = norm(vectorizeData(data),inf);
    data_2norm = norm(vectorizeData(data),2);
    
    disp(['One-norm of the datacube: ' num2str(data_1norm)])
    disp(['Inf-norm of the datacube: ' num2str(data_infnorm)])
    disp(['Two-norm of the datacube: ' num2str(data_2norm)])

    
    % determine sigma (expected energy mismatch between predicted and recorded data)
    sigma = options.relError * data_2norm;
    
    % Convert the specified wavelet physical length to internal windowing operator index
    wavelet_window_start = nt+options.downshift-(floor(options.q_estLength_negT./dt))+1;
    wavelet_window_end = nt+options.downshift+(floor(options.q_estLength_posT./dt)+1);
    wavelet_length = wavelet_window_end - wavelet_window_start + 1;

    % Convert the specified window length in time to windowing operator index
    data_window_start = floor(options.window_startT ./ dt);
    data_window_end = floor(options.window_endT ./ dt);
    
    % Convert anti-causal global muting time to number of samples
    topmute_end = ceil(options.topmuteT ./ dt);
    
%% Construct necessary operators -----------------------
    
    disp('Preparing for Robust EPSI processing...')
    
    % Check whether we need the serial or distributed version of multiple prediction operator
    if options.parallel
        opEPSI = @opEPSI_dist;
    else
        opEPSI = @opEPSI_serial;
    end
    
    % If data-domain time-windowing during the initial spike matching is specified, construct the operator to do so
    if (data_window_start > 0) && (data_window_end > 0)
        OP_TIMEWIND = opTimeWindow(data, data_window_start, data_window_end);
        if options.parallel
            OP_TIMEWIND = oppKron2Lo(OP_TIMEWIND, opDirac(nr*ns));
        end
    else
        OP_TIMEWIND = opDirac(nt*nr*ns);
    end
    
    % Construct the initial EPSI operator
    OP_EPSI = opEPSI(data,q,topmute_end,dt,options.useOblique);
    
    % Pre-allocate a distributed empty solution if parallel processing is on
    if options.parallel
        preallocatedEmptySolutionVector = distVectorize(distributed.zeros(nr,ns,nt));
    end
    
    % Construct sparsity operator
    if options.useSparsity
        disp('The use of sparsifying transform is ON')
        [OP_SPARSITY, preallocatedEmptySolutionVector] = opSparsity(nt, nr, ns, options);
    else
        disp('The use of sparsifying transform is OFF')
        OP_SPARSITY = opDirac(nt*nr*ns);
    end

%% Begin Main Program =========================================================

if options.invertQ
    %% Joint inversion for G and Q
    
    % Initialize variables for estimates of source signature
    OPPAD_S = opPadTop(q(:),nt_conv);
    q_est = zeros(nt_conv,1);
    q_total = zeros(nt_conv,1);
    q_est(:,1) = OPPAD_S * q(:);
    q_total(:,1) = OPPAD_S * q(:);
    
    % Initialize variables for optimization iterations
    alterloop_itercount = 1; % count the number of loops over both variables
    min_spgiter_between_matching = options.minL1iters; % minimum number of iterations for P before matching for Q
    prev_residual = data_2norm;
    residual = data_2norm;
    residual_log = [];
    
    %% Initial spike matching iteration (getting a sparse approximation of the auto-correlation with spgl1)
    disp(' ')
    disp('CURRENTLY WORKING ON: Initial sparse estimation of primaries from autocorrelation')
    
    % In this mode, source signature q is initialized with zero. 
    % Therefore, OP_EPSI contains only the multiple-prediction term.
    % Rename OP_EPSI to reflect this:
    OP_MULTPRED = OP_EPSI;
    
    % Try to pick out the first reflected event. An accurate estimate here will greatly help ensure the overall quality of the solution
    x = EPSI_SLIM_pickSpikyFirstEvent(data, nt, nr, ns, OP_MULTPRED, OP_TIMEWIND, wavelet_length, options);
    
    % set this initial estimate aside for possible QC if we specify it should be saved
    if options.savepickedspikes
        % warning, will gather distributed solution to the head node
        if options.parallel
            x = distVec2distArray(x, [nr ns nt]);
            x = distSort_timeSlice2shot(x);
            x = undist(x);
        else    
            x = reshape(x, [nt,nr,ns]);
        end
        initial_picked_spikes = x;
        save(options.sol_file,'initial_picked_spikes');
    	disp('Recording initial picked spikes... solution mat file saved successfully')
    	return
    end
    
    clear OP_MULTPRED
    
    
    %% == Main alternating optimization loop ==
    while (total_itercount < options.maxTotalIter) && (alterloop_itercount < options.maxOuterIter)
        
        % make sure x is in physical domain (already is for first iteration after spike picking)
        if alterloop_itercount == 1
            x_phys = x;
        else
            x_phys = OP_SPARSITY' * x;
        end
        
        %% == Wavelet inversion phase ==
        % Check if we actually perform the wavelet matching, judging by sufficient descent in objective
        if alterloop_itercount == 1 % match if first iter
            do_wavelet_match = true;
            wavelet_reg_lambda = options.wavelet_reg_lambda; % L2-energy damping in least-squares inversion for Q
        elseif residual < prev_residual * 0.9 % match if sufficient descent
            do_wavelet_match = true;
            prev_residual = residual;
            wavelet_reg_lambda = options.wavelet_reg_lambda; % L2-energy damping in least-squares inversion for Q
        else
            do_wavelet_match = false;
            disp('Skipping next wavelet matching phase: insufficient descent')
        end
        
        if do_wavelet_match    
            disp(' ')
            disp(['CURRENTLY WORKING ON: Alternating optimization loop ' int2str(alterloop_itercount) ' wavelet inversion phase'])

            % define data for Wavelet inversion, subtract previous estimate from current residual
            b = [vectorizeData(data) - (OP_EPSI * x_phys)];
            clear OP_EPSI
            
            % match for global source signature term
            [update_q update_q_unwindowed] = EPSI_SLIM_invertForQ(b, opEPSI, x_phys, nt, nr, ns, nt_conv, dt, q, wavelet_window_start, wavelet_window_end, wavelet_reg_lambda, options);
            
            % apply update to form new wavelet estimate
            q_total(:,alterloop_itercount+1) = update_q_unwindowed(:); % mostly for record keeping
            q_est(:,alterloop_itercount+1) = update_q(:) + q_est(:,alterloop_itercount);
            
            clear b
            clear x_phys
        else
            % insufficnet descent in residual from P inversion, skip wavelet matching
            q_total(:,alterloop_itercount+1) = zeros(length(q_est(:,alterloop_itercount)),1);
            q_est(:,alterloop_itercount+1) = q_est(:,alterloop_itercount);
        end
        
        
        %% == Primary Green's Function inversion phase via SPGl1 ==
        disp(' ')
        disp(['CURRENTLY WORKING ON: Alternating optimization loop ' int2str(alterloop_itercount) ' primary IR inversion phase'])
        
        % If first iteration, use zeros as initial guess
        % Otherwise, start with result from previous iteration
        if alterloop_itercount == 1
            if exist('preallocatedEmptySolutionVector','var') % emptySolutionVector is generated from the domain of distributed operator
                xprev = preallocatedEmptySolutionVector;
            else
                xprev = [];
            end
            tau = 0;
        else
            xprev = x;
            tau = norm(x,1);
        end
        clear x
                
        % make new EPSI operators with the newsource estimate
        OP_EPSI = opEPSI(data,q_est(:,alterloop_itercount+1),topmute_end,dt,options.useOblique);
        
        % The operator to invert is now
        A = OP_EPSI * OP_SPARSITY';
        b = vectorizeData(data);
        
        % Do the next SPGl1 iteration for inverting P
        if (alterloop_itercount == 1)
            options.decTol = 5e-2;
        elseif (alterloop_itercount < 6)
            options.decTol = 5e-3; % cooling the tolerance of spgl1 root-finding for higher iterations
        else
            options.decTol = 5e-4;
        end
        
        % Perform the L1 optimization:
        %      min  |x|_1  s.t. |Ax - b|_2 < sigma
        % using xprev as warmstart, and tau as initial 1-norm score 
        
        [x, residual, residual_list, iters_performed, sigma_reached] = EPSI_SLIM_invertForG_L1(A, b, sigma, xprev, tau, options);
        
        %% == Prepare for next iteration ==
        % Increment all counters and logs
        total_itercount = total_itercount + iters_performed;
        residual_log = [residual_log; residual_list];
        alterloop_itercount = alterloop_itercount + 1;
        disp(['(current cumulated total L1 gradient-step count: ' int2str(total_itercount) ')'])
        clear A
        clear xprev
        
        % Check exit condiiton
        if (residual <= sigma) || sigma_reached % quit if tolerance level in residual reached
            disp('relative sigma reached, exiting alternating optimization loop')
            break
        end
    end
    % Done alternating inversion for P and Q
    
    %% == Do a final polishing step on the wavelet matching if desired
    % this helps for a final bit of clean up in the accuracy of the first arrivals 
    if options.polishWavelet
        disp(' ')
        disp('CURRENTLY WORKING ON: Entering final wavelet polish/debiasing phase ...');
        
        % bring final IR to physical domain
        x_phys = OP_SPARSITY' * x;
        
        % define data for Wavelet inversion, subtract previous estimate from current residual
        b = [vectorizeData(data) - (OP_EPSI * x_phys)];
        clear OP_EPSI
        
        % match for global source signature term
        update_q = EPSI_SLIM_invertForQ(b, opEPSI, x_phys, nt, nr, ns, nt_conv, dt, q, wavelet_window_start, wavelet_window_end, wavelet_reg_lambda, options);
        
        % apply update to form new wavelet estimate
        q_total(:,end+1) = update_q(:); % mostly for record keeping
        q_est(:,end+1) = update_q(:) + q_est(:,end);
        
        clear b
        clear x_phys
    else
        clear OP_EPSI
    end
    
else 
    
    %% This section gets run if invertQ = 0
    % assumption: Q is known, only invert for G
    
    disp(['CURRENTLY WORKING ON: Primary IR inversion with known source signature'])
    
    % directly copy the input wavelet to the output
    OPPAD_S = opPadTop(q(:),nt_conv);
    q_est = OPPAD_S * q(:);
    q_total = OPPAD_S * q(:);
    
    % use SPGl1 to obtain l1 minimizing inversion
    A = OP_EPSI * OP_SPARSITY';
    b = vectorizeData(data);
    
    options.decTol = 5e-3;
    % Perform the L1 optimization:
    %      min  |x|_1  s.t. |Ax - b|_2 < sigma
    
    [x, residual, residual_list, iters_performed, sigma_reached] = EPSI_SLIM_invertForG_L1_noInvertQ(A, b, sigma, [], 0, options);
    
    total_itercount = iters_performed;
    residual_log = residual_list;

end

%% End Main Program ===========================================================

% save sparsity coefficients if instructed to do so (only for debugging)
if options.savecoefmat
    x_coefs = undist(x);
    save(options.coef_file,'x_coefs');
end

% obtain final primary IR in physical domain
if options.useSparsity
    xphys = real(OP_SPARSITY' * x);
else
    xphys = real(x);
end
clear x

%% Final output and cleanup
EPSI_SLIM_writeOutput( opEPSI, data, xphys, q_est, q_total, nt, nr, ns, dt, wavelet_window_start, wavelet_window_end, total_itercount, residual_log, options )

% cleanup
disp('... Finished.')
if nargin == 0
    rmpath(genpath('.'))
end
