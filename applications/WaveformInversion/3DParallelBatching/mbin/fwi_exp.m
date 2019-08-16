function save_file = fwi_exp(fwi_experiment,save_dir,model_dir,data_dir)
% 
% Run parallel 3D FWI with/without batching
% 
% Usage:
%   fwi_exp(fwi_experiment_num,save_dir,model_dir,data_dir);
% 
% Input:
%   fwi_experiment_num     - which fwi_exp##.m to run
%   save_dir               - directory to save results in
%   model_dir              - model directory 
%   data_dir               - data directory
%  
% Curt Da Silva, 2016
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


%% Load experiment params
    if exist('fwi_experiment','var')
        eval(['fwi_exp' num2str(fwi_experiment) ';']);
    else
        error('Need fwi_experiment number');
    end
    
    v2struct(fwi_params_struct);
    
    disp([modelname ' FWI ']);
    disp(['Loaded experiment ' num2str(fwi_experiment) ' params']);
    
    switch modelname
      case 'overthrust_small'
        overthrust_small_preamble;    
      case 'edam'
        edam_preamble;
      otherwise
        error('Unrecognized model');
    end
    
    save_file = [save_dir 'fwi_results_' modelname '_exp' num2str(fwi_experiment) '.mat'];
    
    %% Model setup
    disp('Loading model');
    
    [v,model] = load_model([model_dir model_file], [nx_s,ny_s,nz_s],modelname,entire_model);
    nx = model.n(1); ny = model.n(2); nz = model.n(3);
    
    ns = initialv_smooth;
    v0 = reshape(opKron(opSmooth(nz,ns),opSmooth(ny,ns),opSmooth(nx,ns))*vec(v),model.n);
    minv = min(vec(v)); maxv = max(vec(v));
    
    %% Set up source/receiver grid 
    % options for solving the helmholtz equations
    lsopts = LinSolveOpts();
    lsopts.tol = 1e-6;
    lsopts.maxit = 10000;
    lsopts.maxinnerit = 5;
    lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
    lsopts.precond = LinSolveOpts.PREC_MLGMRES;
    
    % pdefunc options
    pdeopts = PDEopts();
    pdeopts.helm_dt = model.d;
    pdeopts.debug_mode = false;
    pdeopts.helm_pml_max = 10;
    pdeopts.helm_mat_free = true;
    pdeopts.helm_scheme = PDEopts.HELM3D_OPERTO27;
    pdeopts.numcompsrc = numsrc_mt;
    
    opts = struct;
    opts.pdefunopts = pdeopts;
    opts.lsopts = lsopts;
    opts.subsample_model = coarsen_model;
    opts.disp_progress = false;
    
    model = load_geometry(model,modelname,numsrc);
    
    % frequencies in Hz
    model.freq = frequencies;
    
    nrx = length(model.xrec); nry = length(model.yrec);
    nfreq = length(model.freq);
    nsx = length(model.xsrc); nsy = length(model.ysrc);
    nsrc = nsx*nsy;
    nrec = nrx*nry;
    
    Q = speye(nsrc);
    
    minfunc_opts = struct;
    minfunc_opts.maxIter = nopt_itr-1;
    minfunc_opts.optTol = -1;
    minfunc_opts.obj_aux = false;
    %%
    if nfreq==1
        freq_batch = {1};
    else
        freq_batch = num2cell(1:nfreq,1);
    end
    for i=2:nfreq_cycles
        fbase = num2cell(1:nfreq,1);
        freq_batch = {freq_batch{:},fbase{:}};
    end
    vhist = cell(length(freq_batch),1);
    vest = vec(v0);
    
    for j=1:length(freq_batch)
        fbatch = freq_batch{j};
        % Select only sources at this frequency batch
        srcfreqmask = false(nsrc,nfreq);
        srcfreqmask(:,fbatch) = true;
        opts.srcfreqmask = srcfreqmask;
        disp(['Loading data for freqs ' num2str(model.freq(fbatch))]);    
        freqslice_path = @(fidx) [data_dir modelname '-data-freq-' num2str(model.freq(fbatch(fidx))) '.rsf'];
        D = dload_fslices(freqslice_path,[nrec,nsrc,length(fbatch)]);
        
        switch optimization
          case 'lbfgs_full'
            [obj,vest_sub,comp_grid] = misfit_setup(vest,Q,D,model,opts);
            tic;
            vest_sub = minConf_TMP(obj,vest_sub,vec(minv),vec(maxv),minfunc_opts);
            disp(['Optimization time ' num2str(toc) 's']);
          case {'sto_lbfgs','sto_spg'}            
            opts.batch_mode = true;
            [obj,vest_sub,comp_grid] = misfit_setup(vest,Q,D,model,opts);
            
            frugal_opts = struct;
            spmd, bmax = size(getLocalPart(D),2); end
            frugal_opts.bmax = [bmax{:}];
            frugal_opts.maxIter = nopt_itr;       
            frugal_opts.innerit = innerit;
            switch optimization
              case 'sto_lbfgs'
                frugal_opts.solver = 'lbfgs';
              case 'sto_spg'
                frugal_opts.solver = 'spg';
            end            
            disp('Solving optimization problem');
            tic;
            vest_sub = minfunc_frugal(obj,vest_sub,minv,maxv,frugal_opts);
            disp(['Optimization time ' num2str(toc) 's']);   
          otherwise
            error('Unrecognized optimization method');
        end
        vest = comp_grid.to_fine*vest_sub;
        vhist{j} = vest;
        save(save_file,'v0','vhist','vest');
    end    
end