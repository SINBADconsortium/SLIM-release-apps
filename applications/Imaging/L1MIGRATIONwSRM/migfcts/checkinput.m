function options = checkinput(options)
% Syntax:
% options = checkinput(options)
%
% Note:
% This function is exclusively called by Migration_with_SRM. The only purpose to
% separate it is to make the main driver more human-readable.
%
% Description:
% This function checks defaults and conflicts in your input arguments.
%
% Input list:
% options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
%
% Output list:
% options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
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

disp('===Now start checking input arguments.')

%% check input
% input data and source files will be checked when they are read
% by "read_MigData"

% model
if isempty(options.vel_bg)
    error('Fatal: no file containing the background velocity model is present.')
end
% unit of velocity
if isempty(options.unit_vel)
    options.unit_vel = 'm/s';
    disp('Default velocity unit: m/s. Change if this is wrong.')
end

%% default geometry
% receiver X position
if isempty(options.x_rcv_grid)
    error('Fatal: lateral receiver positions not defined.')
end
% receiver Z position
if isempty(options.z_rcv_grid)
    error('Fatal: receiver depth not defined.')
end
% source grid X position
if isempty(options.x_src_grid)
    options.x_src_grid = options.x_rcv_grid;
    disp('Default: source positions are the same as receivers.')
end
% source grid Z position
if isempty(options.z_src_grid)
    options.z_src_grid = options.z_rcv_grid;
    disp('Default: source depth is the same as receivers.')
end
% source location
if isempty(options.x_src)
    options.x_src = options.x_src_grid;
    disp('Default: source presents at every source grid points.')
end

%% check (de)migration parameters if not in "test_readdata" mode
if ~options.test_readdata
    
    %% frequency domain parameters
    % upper frequency boundary
    if isempty(options.freq_upper)
        error('Fatal:Upper boundary of frequency undefined.')
    end
    % lower frequency boundary
    if isempty(options.freq_lower)
        disp('Lower boundary of frequency undefined. Use 0Hz.')
        options.freq_lower = 0;
    end
    
    %% default wavelet type and source polarity
    if isempty(options.wavelet_type)
        options.wavelet_type = 'IP';
    end
    if isempty(options.src_polarity)
        options.src_polarity = 'di';
    end
    
    %% check for "test_model" mode
    if options.test_model
        if isempty(options.source_type)
            % by default type=1, use ricker wavelet
            options.source_type = 1;
            options.data_type = 1;
            options.wavelet_type = 'RK';
        end
        options.withShotRenewal = 0;
	options.withFreqRenewal = 0;
    end
    
    %% check for "model_data" mode
    if options.model_data
        if isempty(options.model_type)
            error('Model type needs to be designated.')
        end
        options.withShotRenewal = 0;
	options.withFreqRenewal = 0;
    end
    
    %% check for "make_linear_data" mode
    if options.make_linear_data
        options.withShotRenewal = 0;
	options.withFreqRenewal = 0;
    end
    
    %% check for migration mode
    % RTM
    if options.runRTM
        options.withShotRenewal = 0;
	options.withFreqRenewal = 0;
	options.sparse_base = 'dirac';
        if isempty(options.RTM_file_name)
            options.RTM_file_name = 'RTM_results.mat';
        end
    end
    % deconvolutional imaging
    if options.runDeconv
        options.withShotRenewal = 0;
	options.withFreqRenewal = 0;
	options.sparse_base = 'dirac';
        if isempty(options.smooth_diameter)
            error('Fatal: smooth diameter needs to be defined for deconvolutional imaging.')
        end
        if isempty(options.deconv_file_name)
            options.deconv_file_name = 'deconv_results.mat';
        end
    end
    % in inversion: use deconvolutional preconditioner
    if options.use_decon_precon
        if isempty(options.smooth_diameter)
            error('Fatal: smooth diameter needs to be defined for deconvolutional preconditioning.')
        end
    end
    % in inversion: if use L2 solver, use dirac basis
    if options.L2_solver
        options.sparse_base = 'dirac';
    end
    
    %% renewal
    % frequency renewal
    if ~options.withFreqRenewal
        if isempty(options.frequency) && isempty(options.nb_freq)
            disp('Frequency not defined. Use all frequencies.')
        end
    else
        if isempty(options.nb_freq)
            error('Fatal: number of frequencies to use is required if using frequency renewal.')
        end
        options.frequency = [];
    end
    % source renewal
    if ~options.withShotRenewal
        if isempty(options.source_encoding_mat) && isempty(options.nb_super_shot)
            disp('Source encoding term is not defined. Use all sequential sources.')
        end
    else
        if isempty(options.nb_super_shot)
            error('Fatal: if using shot renewal, number of supershots needs to be defined.')
        end
        options.source_encoding_mat = [];
    end
    
    %% default output
    if isempty(options.results_folder)
        error('Fatal: results folder not defined.')
    elseif ~strcmp(options.results_folder(end),'/')
        options.results_folder = [options.results_folder,'/'];
    end
    
    if isempty(options.preview_folder)
        options.preview_folder = options.results_folder;
    elseif ~strcmp(options.preview_folder(end),'/')
        options.preview_folder = [options.preview_folder,'/'];
    end
    
    if options.save_evolve && isempty(options.evolve_name)
        options.evolve_name = 'model_updates_file.mat';
    end
    
    if options.save_solmat && isempty(options.solmat_name)
        options.solmat_name = 'solution_file.mat';
    end
    
    if isempty(options.preview_name)
        options.preview_name = 'preview_file.mat';
    end
end

disp('Everything looks fine so far.')