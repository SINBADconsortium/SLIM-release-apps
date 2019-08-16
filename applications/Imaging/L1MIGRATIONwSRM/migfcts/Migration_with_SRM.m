function Migration_with_SRM(varargin)
% Description:
% Sparsity promoting migration with surface related multiples (total upgoing
% wavefield).
%
% Main functionality:
% 	Least-squares two-way migration with sparsity promotion from primaries/surface
% 	related multiples/total upgoing wavefield/EPSI inverted Green's function,
%	depending on the input field "data_type".
%
% Bonus functionalities:
% 	1.) For synthetic examples, test whether the background model is smoothed
%		adequately/whether the Born operator passes the gradient test/whether the
%		absorbing boundary is thick enough/etc. Type "help test_model" in matlab
%		command line for more info. This requires the auxiliary input "vel_true"
%		that you should have for a synthetic example.
%	2.) For synthetic examples, make linearized data. What kind of linearized
%		data you want to make (Green's function/primary/multiple/total data)
%		depends on the input field "data_type". Type "help make_linear_data" in
%		matlab command line for more info. This requires the auxiliary input
%		"model_pert".
%	3.) Cross-correlation reverse time migration. From either primaries/multiples/
%		total data, depending on the input field "data_type".
%	4.) Least-squares two-way migration WITHOUT sparsity promotion. LSQR is used
%		as solver. Damped least-squares algorithm is not used since choosing the
%		right damping factor is non-trivial, which involves a lot of parameter
%		tuning and tweaking. Note that no sparsity transform is used since no
%		sparse regularization is used. Set the input field "L2_solver" to 1 if
%		you want to see how L2 works.
%
% Input: variable length input argument list. For a detailed list, type "edit
% Migration_with_SRM" in matlab command line. You can then see the full list of
% input variables and their default values. Note that a few fields have default
% values while most are empty.
%
% Output: results will be saved to file.
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

options = process_options_to_struct(varargin, ...
    'vel_bg',               '',...          % input: background velocity, variable name has to be 'vel_bg' if you use mat files
    ...                                     % unit: velocity unit, can be m/s or km/s, which needs to be given in field "unit_vel"
    'vel_true',             '',...          % auxiliary input: true velocity model
    ...                                     % unit: velocity unit, can be m/s or km/s, which needs to be given in field "unit_vel"
    ...                                     % NEEDED for testing vel_bg/ FD modeling operator/Born modeling operator
    'model_pert',           '',...          % auxiliaty input: model perturbation
    ...                                     % unit: in slowness square: 1/(unit_vel)^2
    ...                                     % NEEDED for making linearized data
    'unit_vel',             '',...          % indicates whether 'km/s' or 'm/s' is used
    'constant_density',     1,...           % curretly only support constant-density modelling
    'make_linear_data',     0,...           % make linearized data, the data will be:
    ...                                     % 1. Green's function
    ...                                     % 2. primaries
    ...                                     % 3. multiples
    ...                                     % 4. total data
    ...                                     % depending on the input field "data_type" accordingly.
    ...                                     % The program will quit after finishing making linearized data.
    ...                                     % Type "help make_linear_data" in matlab command line for more info.
    'linear_data_name',     '',...          % name of linearized data, if empty, use 'linear_data.mat'
    'length_boundary',      300,...         % length of the absorbing boundary, in meters.
    ...                                     % Default value is 300m, which corresponds to the wavelength of 5Hz sinusoid in water.
    'test_model',           0,...           % test mode for the smooth and true models:
    ...                                     % 1. whether the background velocity/density models are smooth enough not to give rise to any reflection
    ...                                     % 2. whether the absorbing boundary is wide enough for wave at boundary to die out
    ...                                     % 3. difference between the linearized data and the nonlinear residual [F(m)-F(m0)]
    ...                                     % The program will quit after tests finish.
    ...                                     % Type "help test_model" in matlab command line for more info.
    'test_file_name',       '',...          % file name of model-testing mode, if empty, use 'test_model_data.mat'
    'model_data',           0,...           % make background data mode
    'model_type',           [],...          % 1: vel_bg; 2: vel_true
    'model_direct_wave',    1,...           % whether to model direct wave
    'model_data_name',      '',...          % file name of the modelling data, if empty, use "model_data.mat"   
    'test_readdata',        0,...           % if set to 1, only test whether data is correctly read
    'runRTM',               0,...           % only run cross-correlation migration (a single RTM)
    ...                                     % The program will quit after finishing running a RTM.
    'RTM_file_name',        '',...          % file name to save the cross-correlation migration result
    ...                                     % if empty, will be 'RTM_results'
    'runDeconv',            0,...           % only run deconvolutional migration
    ...                                     % The program will quit after finishing running deconvolutional migration
    'deconv_file_name',     '',...          % file name to save the deconvolutional image
    ...                                     % if empty, will be "deconv_results"
    'use_decon_precon',     0,...           % use deconvolutional imaging operator as a preconditioner
    'smooth_diameter',      [],...          % smooth diameter in deconvolutional imaging operator, in meters
    ...                                     % when < 1, it is the Wiener filter parameter
    'L2_solver',            0,...           % using LSQR as solver for inversion. To solve the problem:
    ...                                     % minimize ||Kdm-data||_2
    ...                                     % where K is the Born operator
    ...                                     % By default, SPGl1 is used as solver, which is to solve:
    ...                                     % minimize ||x||_1, subject to ||KCx-data||_2 < sigma
    ...                                     % where "sigma" is the L2 tolarance for data misfit
    ...                                     % Damped least-squares is NOT used since choosing the damp factor is non-trivial.
    'L2_renewal_iter',      10,...          % if use renewal for L2 solver, how frequently (after how many iterations)
    ...                                     % do we draw a new subsampling operator
    'L2_projection',        0,...           % using L2 projection instead of L1 
    'robust_misfit',        0,...           % using robust misfit:
    ...                                     % 0: SPGl1 default: L2 norm
    ...                                     % 1: huber
    ...                                     % 2: student T
    'robust_param',         [],...          % either huber parameter or nu for student T
    'EPSI_solmat',          '',...          % full path of the solution file of EPSI
    ...                                     % This solution file contains:
    ...                                     % 1. EPSI inverted Green's function
    ...                                     % 2. EPSI inverted primaries
    ...                                     % 3. EPSI inverted multiples
    ...                                     % 4. EPSI inverted source wavelet
    ...                                     % 5. time-padded (pad # of time samples to the power of 2) total data
    ...                                     % NOTE: This input has lower priority than input field "separate_input".
    ...                                     % If "separate_input" field is 1, then those data will be read.
    'separate_input',       0,...           % whether to read data individual input files
    ...                                     % If set to 1, input field "EPSI_solmat" will be ignored.
    ...                                     % Instead, the following six fields will be needed
    'input_data',           '',...          % file name of input data
    'data_variable_name',   '',...          % variable name of input_data if is matlab file
    'input_point_source',   '',...          % file name of point source Q
    'Q_variable_name',      '',...          % variable name of point source if is matlab file
    'input_areal_source',   '',...          % file name of areal source P
    'P_variable_name',      '',...          % variable name of areal source if is matlab file
    'input_endian',         [],...          % endian for input "SU" or "binary" files
    'data_type',            [],...          % type of data you want to migrate from
    ...                                     % 1: primaryIR/Green's function
    ...                                     % 2. primaries
    ...                                     % 3. multiples
    ...                                     % 4. total data(primaries+multiples)
    'source_type',          [],...          % type of sources in the Born operator
    ...                                     % if empty, will be the same as "data_type"
    ...                                     % 1. impulse source (to map Green's function)
    ...                                     % 2. wavelet estimated with EPSI: Q (to map primaries)
    ...                                     % 3. areal source: total data -P (to map the multiples)
    ...                                     % 4. Q-P (to map the total data)
    'wavelet_type',         '',...          % type of wavelet to simulate the impulse response, can be:
    ...                                     % 1. a real impulse source by setting to 'IP'
    ...                                     % 2. a band-limited Ricker wavelet by setting to 'RK'
    ...                                     % NOTE: even for real impulse source, a Tukey window with is applied.
    'TukeyWin_alpha',       0.4,...         % if wavelet type is IP, the tapering parameter of the Tukey window.
    'taper_low_freq',       0,...           % if wavelet type is IP, whether the Tukey window tapers both low
    ...                                     % and high frequencies or only the high frequency
    'src_polarity',         '',...          % polarity of source to simulate the impulse response, can be:
    ...                                     % 1. dipole source ('di') to include source ghost
    ...                                     % 2. use monople source ('mo') not to
    'wavelet_shift',        0,...           % shift of the wavelet, in seconds
    'use_oblique',          0,...           % whether to correct for 2D effect, should be kept the same as in EPSI
    'scale_coeff',          1,...           % scale data to a proper amplitude for higher numerical precision
    'model_pert_out',       [],...          % optional output: model perturbation
    ...                                     % if not empty, a result file will be written onto disk
    ...                                     % file format depends on the suffix of this input field
    ...                                     % can be SU/RSF/MAT/binary files, given ".su/.rsf/.mat/.bin files"
    'results_folder',       [],...          % directory to put output files
    'preview_folder',       [],...          % directory to put preview files
    ...                                     % optional. It will be the same as "results_folder" if is emtpy.
    'save_evolve',          1,...           % whether to save solution updates
    'evolve_name',          [],...          % name of file to store the updates
    ...                                     % can only be matlab files
    ...                                     % is 'model_updates_file.mat' by default
    'save_solmat',          1,...           % whether to save solution and diagnostic information (solmat)
    'solmat_name',          [],...          % name of file to store solmat
    ...                                     % can only be matlab files
    ...                                     % is 'solution_file.mat' by default
    'preview_name',         [],...          % name of the preview file
    ...                                     % output of preview file is MANDATORY since all other outputs are optional
    ...                                     % the preview file contains
    ...                                     % 1. model perturbation "dm"
    ...                                     % 2. a single RTM result "dm_RTM"
    ...                                     % 3. diagnostic info: input arguments "options"
    ...                                     % 4. diagnostic info: modelling parameters "model_para"
    ...                                     % is 'preview_file.mat' by default
    'show_result',          1,...           % whether to display preview after the program finishes
    ...                                     % only available if inversion is done.
    'color_data',           [],...          % colormap scheme when data is shown: gray by default
    'color_model',          [],...          % colormap scheme when model is shown: gray by default
    'iter',                 100,...         % number of iterations
    'dt',                   [],...          % sampling interval in time, in seconds
    'nt',                   [],...          % number of samples in time
    'dx',                   [],...          % horizontal discretization of model, in meters
    'dz',                   [],...          % vertical discretization of model, in meters
    'nx',                   [],...          % number of horizontal grid points of model
    'nz',                   [],...          % number of vertical grid points of model
    'x_rcv_grid',           [],...          % horizontal-coordinates of the receiver grid points, in meters
    'z_rcv_grid',           [],...          % depth of the receiver grid points, in meters
    'x_src_grid',           [],...          % horizontal-coordinates of the source grid points, in meters
    'z_src_grid',           [],...          % depth of the source grid points, in meters
    'x_src',                [],...          % horizontal-coordinates of the sources, in meters
    ...                                     % note this is different from that of source grid points
    ...                                     % sources can be absent on some grid points
    'freq_lower',           [],...          % lower frequency boundary
    'freq_upper',           [],...          % upper frequency boundary
    'freq_peak',            [],...          % peak frequency
    'mute_model',           0,...           % whether to mute the source layer from the model
    'min_waterdepth',       [],...          % depth of water layer,
    ...                                     % for muting sources in image space
    'withFreqRenewal',      0,...           % quit Pareto curve upon reaching it, then draw a new bunch of frequencies
    ...                                     % 1. if is 1: only "nb_freq" is needed, "frequency" is ignored
    ...                                     % 2. if is 0 then: "frequency" is used if not empty, "nb_freq" is ignored;
    ...                                     %	if "frequency" is empty but "nb_freq" is not, "nb_freq" is used to choose nb_freq random frequencies
    ...                                     %	if both "frequency" and "nb_freq" are empty, use all frequency
    'withShotRenewal',      0,...           % quit Pareto curve upon reaching it, then draw a new bunch of shots
    ...                                     % 1. if is 1: only "nb_super_shot" is used, "source_encoding_mat" is ignored;
    ...                                     % 2. if is 0 then: "source_encoding_mat" is used if not empty,"nb_super_shot" is ignored
    ...                                     %	if "source_encoding_mat" is emtpy but "nb_super_shot" is not, "nb_super_shot" is used to make a random Gaussian matrix with "nb_super_shot" columns
    ...                                     %	if both "source_encoding_mat" and "nb_super_shot" are empty, use all sequential shots
    'nb_freq',              [],...          % number of frequencies
    'frequency',            [],...          % frequencies to use if not using frequency renewal
    'nb_super_shot',        [],...          % number of supershots
    'encoding_mat_type',    'gaussian',...  % type of encoding matrix: can be "gaussian" or "dirac"
    'source_encoding_mat',  [],...          % source encoding matrix to use if not using shot renewal
    'initial_guess',        [],...          % initial guess for the solution
    'sparse_base',          'curvelet',...  % can be set to be either 'curvelet' (recommended) or 'wavelet(daubechies/haar)' or 'dirac'
    ...                                     % the sparsifying domain where model perturbation is sparse
    'curvelet_type',        'WRAP',...      % type of curvelet boundary condition, can be 'WRAP' or 'ME'
    'rel_error',            0,...           % relative error, 0 by default, meaning solving a BP problem
    'real_flag',            1,...           % apply a real restriction if the output is real, in this case yes.
    'depth_precon',         1 ...           % use depth preconditioner
    );

%% check input: defaults and conflicts
options = checkinput(options);

%% load data, models and reshape to correct size
[data, source, file_dt, model_bg, model_true, model_pert] = read_MigData(options);

% determine dt
if ~isempty(options.dt) && (options.dt ~= 0)
    dt = options.dt;
elseif ~isempty(file_dt) && (file_dt ~= 0)
    dt = file_dt;
else
    error('dt is never specified, neither in the data file nor in the options')
end
options.dt = dt;

% update other related information
[options.nz options.nx] = size(model_bg.vel);
if ~isempty(data)
    options.nt = size(data,1);
elseif ~isempty(source.P)
    options.nt = size(source.P,1);
end

%% show plots if in read data mode
if options.test_readdata
    disp('===Reading data mode: displaying data/model only.')
    figure;
    plot_utils.show_data(data(:,:,floor(size(data,3)/2)),options);
    title('Preview of data: one shot gather.')
    figure;
    plot_utils.show_model(model_bg.vel,options);
    title('Preview of background model: velocity.')
    if ~isempty(model_true.vel)
        figure
        plot_utils.show_model(model_true.vel,options);
        title('Preview of true model: velocity.')
    end
    if ~isempty(model_pert.slowness2)
        figure
        plot_utils.show_model(model_pert.slowness2,options);
        title('Preview of model perturbation: slowness square.')
    end
    if ~isempty(source.Q)
        figure
        plot_utils.show_data(source.Q, options);
        title('Preview of source wavelet.')
    end
    if ~isempty(source.P)
        figure
        plot_utils.show_data(source.P(:,:,floor(size(source.P,3)/2)),options);
        title('Preview of areal source.')
    end
    disp('===Finish reading data. Quit now.')
    return
end

disp('===Processing acquisition geometry/frequencies/sources/etc.')

%% geometries
dx = options.dx;                        % horizontal grid distance
dz = options.dz;                        % vertical grid distance
nx = options.nx;                        % number of samples of the model horizontally
nz = options.nz;                        % number of samples of the model vertically
zs_grid = options.z_src_grid;           % z-coordinates of source grid points
zr_grid = options.z_rcv_grid;           % z-coordinates of receivers grid points
xd = (nx-1)*dx;                         % Real distance between horizon grid points
zd = (nz-1)*dz;                         % Real distance between depth grid points
xr_grid = options.x_rcv_grid;           % x-coordinates of receiver grid points
xs_grid = options.x_src_grid;           % x-coordinates of source grid points
xs = options.x_src;                     % x-coordinates of sources
nrec = length(xr_grid);                 % number of receivers
nshot = length(xs);                     % number of sources
ns_grid = length(xs_grid);              % number of source grid points
nbx = floor(options.length_boundary/dx);
nbz = floor(options.length_boundary/dz);

% mute model parameter
mute_model = options.mute_model;
min_waterdepth = options.min_waterdepth;
if mute_model~=0
    if isempty(min_waterdepth) || min_waterdepth==0
        mute_model = 0;
    else
        mute_model = 1;
        topmute_model = floor(min_waterdepth/dz);
    end
end

% processing for dipole sources: pad model and sink source/receivers
src_polarity = options.src_polarity;
if strcmp(src_polarity,'di')
    % pad background model
    padsize = ceil(zs_grid/dz)+1;
    vel_pad = repmat(model_bg.vel(1,:),[padsize,1]);
    model_bg.vel = [vel_pad;model_bg.vel];
    % pad model perturbation if known
    if ~isempty(model_pert.slowness2)
        pert_pad = repmat(model_pert.slowness2(1,:),[padsize,1]);
        model_pert.slowness2 = [pert_pad;model_pert.slowness2];
    end
    % pad true model if known
    if ~isempty(model_true.vel)
        vel_pad_true  = repmat(model_true.vel(1,:),[padsize,1]);
        model_true.vel = [vel_pad_true;model_true.vel];
    end
    % depth dimension of model
    nz = nz + padsize;
    % depth offset of model
    oz = -dz*padsize;
    % dipole source
    zs_grid = [-zs_grid zs_grid];
    % mute depth
    if mute_model
        topmute_model = topmute_model + padsize;
    end
else
    padsize = 0;
    oz = 0;
end

%% processing time-frequency relationships
nt = options.nt;
freq_lower = options.freq_lower;
freq_upper = options.freq_upper;

% the Nyquist frequency
if mod(nt,2) == 1
    nf_nyquist = ceil(nt/2);
    has_nyquist = 0;
else
    nf_nyquist = nt/2+1;
    has_nyquist = 1;
end
df = 1/(nt*dt);
frequency_nyquist = 0:df:df*(nf_nyquist-1);
frequency_nyquist = vec(frequency_nyquist);
if freq_upper > frequency_nyquist(end)
    disp(['Upper frequency boundary LARGER than Nyquist. Use Nyquist ' ...
          'frequency instead.'])
    freq_upper = frequency_nyquist(end);
end

% frequency_full is the full range of frequencies in [freq_lower, freq_upper]
frequency_full = df*ceil(freq_lower/df):df:df*floor(freq_upper/df);
frequency_full = vec(frequency_full);
if frequency_full(1) == 0
    % does not contain 0 frequency
    frequency_full = frequency_full(2:end);
end
nf_full = size(frequency_full,1);
% related variables
% frequency_full w.r.t Nyquist
fmask_full = get_subset_mask(frequency_nyquist,frequency_full);
% frequency_full w.r.t itself
fidx_full = get_subset_mask(frequency_full,frequency_full);
% make Tukey window for IMPULSIVE source
if ~strcmp(options.wavelet_type,'RK')
    taper_sides.left = options.taper_low_freq;
    taper_sides.right = 1;
    opFreqTaper_full = opTukeyWinMask(nf_full,options.TukeyWin_alpha,fidx_full,taper_sides);
else
    opFreqTaper_full = opDirac(nf_full);
end

% choose frequency subset within frequency_full
if options.withFreqRenewal == 0
    % no frequency renewal
    frequency = options.frequency;
    if isempty(frequency)
        % choose frequency if not designated
        if isempty(options.nb_freq)
            % use all frequencies
            frequency = frequency_full;
        else
            % use frequency subset
            if options.nb_freq > nf_full
                disp(['Warning: designated number of frequencies ' ...
                      'more than total number of frequencies. Use ' ...
                      'all frequencies instead.'])
                options.nb_freq = nf_full;
            end
            frequency = choose_freq_subset(frequency_full, options.nb_freq);
        end
    else
        % put designated frequencies onto grid points
        disp('Map given frequencies to grids.')
        frequency = df.*round(frequency./df);
        % check frequency range
        if (frequency(1) < frequency_full(1))
            frequency(1) = frequency(1) + df;
            if (frequency(1) == frequency(2))
                frequency = frequency(2:end);
                disp('Redundant frequency removed.')
            end
        end
        if (frequency(end) > frequency_full(end))
            frequency(end) = frequency(end) - df;
            if (frequency(end) == frequency(end-1))
                frequency = frequency(1:end-1);
                disp('Redundant frequency removed.')
            end
        end
        if (frequency(1) < frequency_full(1)) || (frequency(end) > frequency_full(end))
            error('Input frequencies not in the designated frequency band.')
        end
    end
    frequency = vec(frequency);
    nf = size(frequency,1);
    fprintf('Using %d frequencies. Listed below: \n', nf);
    for freq_count = 1:nf
        fprintf('Frequency (%d): %.4f \n',freq_count, frequency(freq_count));
    end
    % frequency mask w.r.t the frequencies up to nyquist
    fmask = get_subset_mask(frequency_nyquist, frequency);
    % frequency mask w.r.t frequency_full
    fidx = get_subset_mask(frequency_full,frequency);
else
    nf = options.nb_freq;
    fprintf('Using %d frequencies. Frequencies renew with iterations.\n',nf);
    % other parameters will be defined inside migration
end

%% process source info
% source-frequency analogue:
% freq: nyquist-->full-->subset
% source: src_grid-->source-->encoded sources

% find the index where source presents at source grids
src_idx = get_subset_mask(xs_grid,xs);
source_encoding_mat_full = eye(nshot);

% source encoding
if options.withShotRenewal == 0
    source_encoding_mat = options.source_encoding_mat;
    if isempty(source_encoding_mat)
        % make encoding mat if not designated
        if isempty(options.nb_super_shot)
            % use all sequential sources
            source_encoding_mat = source_encoding_mat_full;
        else
            % make an encoding matrix
            if options.nb_super_shot > nshot
                error(['Number of encoded sources should not be larger than number of sequential sources.'])
            end
            source_encoding_mat = make_encoding_mat(nshot,options.nb_super_shot,options.encoding_mat_type);
        end
    else
        % check existing encoding matrix
        if not(isequal(size(source_encoding_mat,1), nshot))
            error('Source encoding matrix has wrong dimension.');
        end
        if size(source_encoding_mat,2) > nshot
            error(['Number of encoded sources should not be larger than number of sequential sources.']);
        end
    end
    nss = size(source_encoding_mat,2);
    fprintf('Using %d shots/supershots.\n',nss);
else
    nss = options.nb_super_shot;
    fprintf('Using %d shots/supershots.\n',nss);
    % other parameters will be defined inside migration
end

%% prepare parameter for impluse response. Note that:
% 1. when peak frequency = 0, an impulse source is used
% 2. when peak frequency ~= 0, a Ricker wavelet is used, which yields a band
% 	limited response.
freq_peak = options.freq_peak;

if strcmp(options.wavelet_type,'RK')   	% Ricker wavelet (for band-limited IR)
    if isempty(freq_peak)
        freq_peak = freq_upper/2.4;     % empirical value for ricker wavelet
        fprintf('Ricker wavelet used. Peak freq by default: %f (%f/2.4)\n',freq_peak,freq_upper);
    else
        fprintf('Ricker wavelet used. Peak freq: %f\n',freq_peak);
    end
else
    freq_peak = 0;                      % impulse wavelet (for flat-spectrum)
    disp('Implusive wavelet used for Greens function.')
end

%% transform data to frequency domain, distribute and vectorize
existP_flag = not(isempty(source.P));
existQ_flag = not(isempty(source.Q));
existI_flag = not(isempty(source.I));
existData_flag = not(isempty(data));

% distribute source.P for memory efficiency
if existP_flag
    % check physics for migrating multiples
    if not(isequal(xs_grid, xr_grid))
        error(['Fatal: sources and receivers must be collocated when  ' ...
               'using total data as source to migrate multiples. ' ...
               'Physics are not correct.'])
    end
    % BEGIN: processing source.P in the same way as is done in EPSI.
    % The following part can be further optimized, however,
    % it is kept this way to be consistent with what is done in
    % EPSI. Otherwise the program can break.
    if options.use_oblique
        opDoObliq = opObliq(dt,nf_nyquist,nrec*nshot);
    else
        opDoObliq = opDirac(nf_nyquist*nrec*nshot);
    end
    % Time domain windowing
    time_taper_frac = 0.1; % ranges as a fraction from 0 to 1
    window_start = 1;
    window_end = nt;
    window_length = window_end-window_start+1;
    taper_length = floor((time_taper_frac/2)*window_length);
    opWindTime = opTimeWindow_Tapered(zeros(nt,1),window_start,window_end,taper_length);
    % Receiver coordiante windowing
    recv_taper_frac = 0.05; % ranges as a fraction from 0 to 1
    window_start = 1;
    window_end = nrec;
    window_length = window_end-window_start+1;
    taper_length = floor((recv_taper_frac/2)*window_length);
    opWindRec = opTimeWindow_Tapered(zeros(nrec,1),window_start,window_end,taper_length);
    % Source coordinate remains untouched
    opWindSrc = opDirac(nshot);
    % Window for datacube: source.P
    opWindCube = opKron(opWindSrc,opWindRec,opWindTime);
    % Fourier transform and apply obliquity factor
    opFT = opKron(opDirac(nrec*nshot),opFFTsym_conv_mask(options.nt));
    opIFT = opFT';
    % apply the windowing and obliquity operators
    source.P = opIFT*opDoObliq*opFT*opWindCube*vec(source.P);
    source.P = reshape(source.P, nt, nrec, nshot);
    % END: processing source.P in the same way as is done in EPSI.
    % shift dimension and vectorize
    source.P = shiftdim(source.P,1);
    source.P = vec(source.P);
    % compute spectrum and apply the taper for impulsive source
    opFT_P = opKron(opFreqTaper_full*opFFTsym_conv_mask(options.nt,fmask_full),opDirac(nrec*nshot));
    sourceP = opFT_P*source.P;
    sourceP = reshape(sourceP, nrec, nshot, nf_full);
    % distribute over frequency
    spmd
        sourceP = codistributed(sourceP, codistributor1d(3));
    end
    % vectorize
    source.P = vec(sourceP);
    % garbage collection
    clear opFT_P sourceP
end

% distribute source.Q for consistency in the code, also for potential
% extension in the future to allow source directivity
if existQ_flag
    % compute spectrum and apply the taper for impulsive source
    opFT_Q = opFreqTaper_full*opFFTsym_conv_mask(options.nt, fmask_full);
    sourceQ = opFT_Q*source.Q;
    % expand Q to be a cube
    temp_Q = eye(ns_grid);
    temp_Q = temp_Q(:,src_idx);
    sourceQ = opKron(opDiag(sourceQ), opDirac(ns_grid*nshot))* ...
        vec(repmat(temp_Q, [1, 1, nf_full]));
    sourceQ = reshape(sourceQ, ns_grid, nshot, nf_full);
    % distribute over frequency
    spmd
        sourceQ = codistributed(sourceQ, codistributor1d(3));
    end
    % vectorize
    source.Q = vec(sourceQ);
    % garbage collection
    clear opFT_Q temp_Q sourceQ
end

% distribute source.I for consistency in the codealso for potential
% extension in the future to allow source directivity
if existI_flag
    % compute spectrum for impulsive source
    sourceI = ones(nf_full,1);
    % expand I to be a cube
    temp_I = eye(ns_grid);
    temp_I = temp_I(:,src_idx);
    sourceI = opKron(opDiag(sourceI), opDirac(ns_grid*nshot))* ...
        vec(repmat(temp_I, [1, 1, nf_full]));
    sourceI = reshape(sourceI, ns_grid, nshot, nf_full);
    % distribute over frequency
    spmd
        sourceI = codistributed(sourceI, codistributor1d(3));
    end
    % vectorize
    source.I = vec(sourceI);
    % garbage collection
    clear temp_I sourceI
end

% distribute data for memory efficiency
if existData_flag
    switch options.data_type
      case 1
        disp('Data type: Greens function.')
      case 2
        disp('Data type: primaries.')
      case 3
        disp('Data type: multiples.')
      case 4
        disp('Data type: total data.')
    end
    % shift dimension and vectorize
    data = shiftdim(data,1);
    data = vec(data);
    % estimate numerical range for optimal precision
    scale_coeff = options.scale_coeff/max(abs(data));
    % compute spectrum
    opFT_data = opKron(opFFTsym_conv_mask(options.nt,fmask_full),opDirac(nrec*nshot));
    data = opFT_data*data;
    data = reshape(data, nrec, nshot, nf_full);
    % distribute over frequency
    spmd
        data = codistributed(data,codistributor1d(3));
    end
    % vectorize
    data = vec(data);
    % garbage collection
    clear opFT_data
end

%% prepare other frequency independent modeling parameters
model_para.o = [oz 0 0];
model_para.d = [dz, dx, 1];
model_para.n = [nz, nx, 1];
model_para.nb = [nbz, nbx, 0];
model_para.zsrc = zs_grid;
model_para.zrec = zr_grid;
model_para.xsrc = xs_grid;
model_para.xrec = xr_grid;
model_para.f0 = freq_peak;
model_para.t0 = options.wavelet_shift;
model_para.xsrc_is = xs;
model_para.nss = nss;

%% test mode: test model
if options.test_model
    disp('===Test model===')
    model_para.freq = frequency;
    model_para.fmask = fmask;
    model_para.fidx = fidx;
    model_para.source_encoding_mat = source_encoding_mat;
    test_model(model_bg.vel, model_true.vel, model_para, source, options);
    disp('Sucessfully finished.')
    disp('===Quit now.')
    return
end

%% modelling data mode
if options.model_data
    disp('===Modelling data===')
    model_para.freq = frequency;
    model_para.fmask = fmask;
    model_para.fidx = fidx;
    model_para.source_encoding_mat = source_encoding_mat;
    if options.model_type == 1
        vel_model = model_bg.vel;
    elseif options.model_type == 2
        vel_model = model_true.vel;
    else
        error('Unknown model type.')
    end
    model_data(vel_model, model_para, source, options);
    disp('Sucessfully finished.')
    disp('===Quit now.')
    return
end

%% make linearized data
if options.make_linear_data
    disp('===Making linearized data mode.')
    model_para.freq = frequency;
    model_para.fmask = fmask;
    model_para.fidx = fidx;
    model_para.source_encoding_mat = source_encoding_mat;
    make_linear_data(model_bg.vel, model_pert.slowness2, model_para, source, options);
    disp('Sucessfully finished.')
    disp('===Quit now.')
    return
end

%% RTM
if options.runRTM
    model_para.freq = frequency;
    model_para.fmask = fmask;
    model_para.fidx = fidx;
    model_para.source_encoding_mat = source_encoding_mat;
    
    disp('Migration mode: reverse time migration.')
    [~, ~, dm_RTM] = est_opMig_norm(model_bg.vel, data, model_para, source, options);
    if strcmp(src_polarity,'di')
        dm_RTM = dm_RTM(padsize+1:end,:);
    end
    writeDataType([options.results_folder, options.RTM_file_name],dm_RTM,'dm_RTM')
    disp('Successfully finished.')
    disp('===Quit.')
    return
end

%% deconv-migration
if options.runDeconv
    model_para.freq = frequency;
    model_para.fmask = fmask;
    model_para.fidx = fidx;
    model_para.source_encoding_mat = source_encoding_mat;
   
    disp('Migration mode: deconvolutional migration.')
    model_para.winsize = [options.smooth_diameter,options.smooth_diameter,1];
    [~, ~, dm_deconv] = est_opPrecon_norm(model_bg.vel, data, model_para, source, options);
    if strcmp(src_polarity,'di')
        dm_deconv = dm_deconv(padsize+1:end,:);
    end
    writeDataType([options.results_folder, options.deconv_file_name],dm_deconv,'dm_deconv')
    disp('Successfully finished.')
    disp('===Quit.')
    return
end

%% inversion
disp('===Start inversion.')

% save evolutional results
opSave = opSaveSnapshot(nz*nx,[options.results_folder, options.evolve_name]);

% model-muting operator
if options.mute_model == 0
    opModelMute = opDirac(nz*nx);
else
    topmute_start = max(ceil(zs_grid(end)/dz)+padsize+1, round(topmute_model/2));
    topmute_end = topmute_model;
    if topmute_start >= topmute_end
        opModelMute = opDirac(nz*nx);
        disp('Water too shallow. Muting is not applied.')
    else
        opModelMute =  opKron(opDirac(nx),opLinMute(nz,topmute_start,topmute_end));
    end
end

% sparsifying synthesis operator. 
if strcmp(options.sparse_base,'wavelet')==1
    C = opKron(opWavelet(nx),opWavelet(nz));
    C = oppad*C';
elseif strcmp(options.sparse_base,'daubechies')==1
    family = 'Daubechies';
    C = opKron(opWavelet(nx,family),opWavelet(nz,family));
    C = oppad*C';
elseif strcmp(options.sparse_base,'haar')==1
    family = 'Haar';
    C = opKron(opWavelet(nx,family),opWavelet(nz,family));
    C = oppad*C';
elseif strcmp(options.sparse_base,'curvelet')==1
    NBscales = max(1,ceil(log2(min(nz,nx)) - 3));
    NBangles = 16;
    Finest = 1;
    Ttype = options.curvelet_type;
    IS_real = 0;
    C = opCurvelet(nz,nx,NBscales,NBangles,Finest,Ttype,IS_real);
    C = C';
elseif strcmp(options.sparse_base,'dirac')==1
    C = opDirac(nz*nx);
else
    error('unknown sparse transform: can be wavelet(daubechies/haar)/curvelet/dirac');
end

% Real Restriction operator
if options.real_flag == 1
    opGetReal = opRealRestriction(nx*nz);
else
    opGetReal = opDirac(nx*nz);
end

% Depth preconditioning
% See this paper (https://slim.gatech.edu/node/6357) for more info
if options.depth_precon == 1
    z = dz:dz:dz*nz;
    opDepth = opKron(opDirac(nx),opDiag(sqrt(z)));
else
    opDepth = opDirac(nx*nz);
end

% set spgl1 parameters
iter = options.iter;
if options.robust_misfit == 1
    misfit_fun = @funHUB;
    params.hub = options.robust_param;
elseif options.robust_misfit == 2
    misfit_fun = @funST;
    params.nu = options.robust_param;
else
    misfit_fun = @funLS;
    params = [];
end

renew_flag = options.withFreqRenewal || options.withShotRenewal;
quitParetoFlag = renew_flag || options.save_evolve;
if options.L2_projection
    opts = spgSetParms('iterations', iter , ...
                       'verbosity', 2 , ...
                       'bpTol', 1e-3 , ...
                       'optTol', 1e-3, ...
                       'quitPareto', quitParetoFlag, ...
                       'decTol', 5e-2, ...
                       'project', @NormL2_project, ...
                       'primal_norm', @NormL2_primal, ...
                       'dual_norm', @NormL2_dual, ...
                       'funPenalty', misfit_fun);
else
    opts = spgSetParms('iterations', iter , ...
                       'verbosity', 2 , ...
                       'bpTol', 1e-3 , ...
                       'optTol', 1e-3, ...
                       'quitPareto', quitParetoFlag, ...
                       'decTol', 5e-2, ...
                       'project', @NormL1_project, ...
                       'primal_norm', @NormL1_primal, ...
                       'dual_norm', @NormL1_dual, ...
                       'funPenalty', misfit_fun);
end

% iteration count, etc
outerloop = 0;
innerloop = 0;
innerloop_list=zeros(iter,1);
residual_list=zeros(2*iter,1);
tau_list=zeros(2*iter,1);
lambda_list=zeros(2*iter,1);
estimated_opMig_norm = 0;
initial_guess_used = 0;
success_flag = 0;

% loop
while (innerloop < iter) && (not(success_flag))
    % frequency renewal
    if options.withFreqRenewal
        model_para.freq = choose_freq_subset(frequency_full, nf);
        model_para.fmask = get_subset_mask(frequency_nyquist, model_para.freq);
        model_para.fidx = get_subset_mask(frequency_full, model_para.freq);
    else
        model_para.freq = frequency;
        model_para.fmask = fmask;
        model_para.fidx = fidx;
    end
    % source renewal
    if options.withShotRenewal
        model_para.source_encoding_mat = make_encoding_mat(nshot, nss, options.encoding_mat_type);
    else
        model_para.source_encoding_mat = source_encoding_mat;
    end
    
    % the scaling is a very rough estimation, so we use the same scaling for different supershots
    if not(estimated_opMig_norm)
        [~, ~, dm_RTM, opMig_norm] = est_opMig_norm(model_bg.vel, data, model_para, source, options);
        scale_opMig = 1/opMig_norm;
        if options.use_decon_precon
            model_para.winsize = [options.smooth_diameter,options.smooth_diameter,1];
            [~, ~, dm_deconv, opPrecon_norm] = est_opPrecon_norm(model_bg.vel, data, model_para, source, options);
        else
            dm_deconv = [];
            opPrecon_norm = 1;
        end
        scale_opPrecon = 1/(opPrecon_norm);
        estimated_opMig_norm = 1;
    end
    % build Born operator and the subsampling operator
    [opMig, opSub] = est_opMig_norm(model_bg.vel, data, model_para, source, options);
    % build the preconditioning operator
    if options.use_decon_precon
        model_para.winsize = [options.smooth_diameter,options.smooth_diameter,1];
        opPrecon = est_opPrecon_norm(model_bg.vel, data, model_para, source, options);
    else
        opPrecon = opDirac(nf*nrec*nss);
    end
    % combine all operators for inversion
    A = scale_opPrecon*(opPrecon*(scale_opMig*(opMig*(opModelMute*(opDepth*(opGetReal*C))))));
    % compute right hand side
    b = scale_opPrecon*(opPrecon*(scale_coeff*(opSub*data)));
    % amplitude correction since we replaced scale_opMig by scale_coeff for higher numerical precision
    options.scale_correction = scale_opMig/scale_coeff;
    % set misfit corresponding to b
    sigma = options.rel_error*norm(b);
    % process warm-start info
    if not(initial_guess_used)
        if length(options.initial_guess) == size(C,2)
            initial_guess = options.initial_guess/options.scale_correction;
            disp('Initial guess used.')
        else
            if ~isempty(options.initial_guess)
                disp('Initial guess NOT used due to dimension mismatch.')
                disp('Possible cause: optimization in a different domain.')
            end
            initial_guess = [];
        end
        tau = norm(initial_guess,1);
        initial_guess_used = 1;
    end
    % inversion with L2 or L1 solver
    if options.L2_solver
        if not(renew_flag)
            options.L2_renewal_iter = options.iter;
        end
        [x, info.flag, info.rNorm, info.iter, info.rNorm2, info.lambda] = ...
            lsqr(A, b, [], options.L2_renewal_iter, [], [], initial_guess);
        outerloop = outerloop + 1;
        innerloop = innerloop + info.iter;
        innerloop_list(outerloop) = info.iter;
        residual_list(innerloop-info.iter+1:innerloop) = info.rNorm2(2:end)*options.scale_correction;
        tau_list(innerloop-info.iter+1:innerloop) = 0;
        lambda_list(innerloop-info.iter+1:innerloop) = info.lambda*options.scale_correction;
        initial_guess = x;
        success_flag = not(info.flag);
    else
        [x,~,~,info]  = spgl1(A,b,tau,sigma,initial_guess,opts,params);
        outerloop = outerloop + 1;
        innerloop = innerloop + info.iter;
        innerloop_list(outerloop) = info.iter;
        residual_list(innerloop-info.iter+1:innerloop) = info.rNorm2*options.scale_correction;
        tau_list(innerloop-info.iter+1:innerloop) = info.xNorm1*options.scale_correction;
        lambda_list(innerloop-info.iter+1:innerloop) = info.lambda;
        initial_guess = x;
        tau = norm(initial_guess, 1);
        infi_loop = isequal(info.iter,0);
        success_flag = (((info.stat >= 1)&&(info.stat <= 4))||(info.stat == 7))||(infi_loop);
    end
    % save evolutional results
    if options.save_evolve == 1
        opSave*opModelMute*opDepth*opGetReal*C*x;
    end
end

% summerize renewal info
renewal_info.innerloop_list=innerloop_list(1:outerloop);
renewal_info.residual_list=residual_list(1:innerloop);
renewal_info.tau_list=tau_list(1:innerloop);
renewal_info.lambda_list=lambda_list(1:innerloop);

% get back physical domain reflectivity
x = options.scale_correction*x;
dm = opModelMute*opDepth*opGetReal*C*x;
dm = reshape(dm,nz,nx);

if strcmp(src_polarity,'di')
    dm = dm(padsize+1:end,:);
    dm_RTM = dm_RTM(padsize+1:end,:);
    if not(isempty(dm_deconv))
        dm_deconv = dm_deconv(padsize+1:end,:);
    end
end

disp('Sucessfully finished.')
disp('===Start writing output.')

%% write output to file
write_MigData(options, dm, dm_RTM, dm_deconv, x, info, renewal_info, model_para);
if options.show_result == 1
    disp_migPreview([options.preview_folder, options.preview_name]);
end
disp('Sucessfully finished.')
disp('===Quit now.')