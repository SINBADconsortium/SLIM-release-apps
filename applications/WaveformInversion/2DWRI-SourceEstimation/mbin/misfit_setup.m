function func = misfit_setup(Q,D,model,params)
% MISFIT_SETUP - Generates a least-squares FWI or WRI objective for use in an optimization algorithm.
% Handles data parallelization along frequencies and sources automatically.
%
% Usage:
%   func = misfit_setup(Q,D,model,params);
%
% Input:
%   Q     - source matrix
%   D     - vectorized, distributed data
%   model - struct with model parameters
%   params
%
%       .wri             - if true (default), computes WRI objective,
%                          if false, computes standard FWI objective
%       .lambda          - tradeoff parameter (scalar) between PDE and data
%                          misfit (for WRI only)
%       .dist_mode       - 'src'     : distributes over sources
%                          'freq'    : distributes over frequencies
%                          'srcfreq' : distributes over sources and frequencies
%
%
%       .freq_index      - if set, uses the specified frequencies (default: 1:length(model.freq)
%       .src_index       - if set, uses the specified sources (default: 1:size(Q,2) )
%       .hessian         - 'sparse' : WRI only, sparse Hessian (default for WRI)
%                          'gn'     : Gauss-Newton Hessian (default for FWI)
%                          'full'   : Full Hessian
%       .obj_per_worker  - if true, each worker has its own objective defined by
%                          the data, source matrix it contains (default: false,
%                          a single objective for all workers, the usual case)
%       .select_src_freq - if true,
%
%
% output:
%   func - objective function used in an optimization procedure
%
%   e.g, evaluating the objective
%     [fk,gk,hk] = func(m);
%   where
%   fk - objective value
%   gk - gradient
%   hk - spot operator or sparse hessian matrix, depending on the value of params.hessian
%
%   If params.obj_per_worker is true, the objectives are evaluated as
%
%   spmd,
%      [fk,gk,hk] = func(m_local);
%   end
%
%   If params.select_src_freq, the objective is evaluated as
%      obj = func(i_src,i_freq);
%      [fk,gk,hk] = obj(m);
%
%
%  Author: Curt Da Silva, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
nsrc = size(Q,2); nrec = length(model.xrec)*length(model.zrec); nfreq = length(model.freq);
nlabs = parpool_size();

assert(nlabs > 0, 'Need available matlab workers');

if exist('params','var')==0, params = struct; end
if ~isfield(params,'wri'), params.wri = true; end
if ~isfield(params,'hessian'), if params.wri, params.hessian = 'sparse'; else params.hessian = 'gn'; end, end
if isfield(params,'dist_mode'), dist_mode = params.dist_mode; else dist_mode = 'freq'; end
if isfield(params,'src_index'), src_index = params.src_index; else src_index = 1:nsrc; end
if isfield(params,'freq_index'), freq_index = params.freq_index; else freq_index = 1:nfreq; end
if isfield(params,'obj_per_worker'), obj_per_worker = params.obj_per_worker; else obj_per_worker = false; end
if isfield(params,'select_src_freq'), select_src_freq = params.select_src_freq; else select_src_freq = false; end

% Restrict sources and frequencies
D = pSPOT.utils.distVec2distArray(D,[nrec,nsrc,nfreq]);
D = pSPOT.utils.distVectorize( D(:,src_index,freq_index) );
model.freq = model.freq(freq_index);
if size(Q,3)==1
    Q = Q(:,src_index);
else
    Q = Q(:,src_index,freq_index);
end
nsrc = length(src_index); nfreq = length(freq_index);
wri_mode = params.wri;

if isfield(params,'FlagUandAlpha')
    FlagUandAlpha = params.FlagUandAlpha;
else
    FlagUandAlpha = 0;
    params.FlagUandAlpha = 0;
end


% Distribute data properly
switch dist_mode
    case 'src'
        D1 = pSPOT.utils.distVec2distArray(D,[nrec,nsrc,nfreq]);
        spmd,
            D1 = redistribute(D1,codistributor1d(3,codistributor1d.unsetPartition,[nrec,nsrc,nfreq]));
        end
        D1 = pSPOT.utils.DistPermute(D1,[1 3 2]) ;
        [D1,loc_size,start_indices,end_indices] = distribute_nd_array(D1,[nrec,nfreq,nsrc],1);
        spmd,
            i_src = start_indices(1):end_indices(1);
            i_freq = 1:nfreq;
        end
        
    case 'freq' % frequencies
        [D1,loc_size,start_indices,end_indices] = distribute_nd_array(D,[nrec,nsrc,nfreq],1);
        spmd,
            i_freq = start_indices(1):end_indices(1);
            i_src = 1:nsrc;
        end
        
    case 'srcfreq' % sources + frequencies
        [D1,loc_size,start_indices,end_indices] = distribute_nd_array(D,[nrec,nsrc,nfreq],2);
        spmd,
            i_src = start_indices(1):end_indices(1);
            i_freq = start_indices(2):end_indices(2);
        end
end

% Set up local model structs, source matrices
spmd,
    model_loc = model;
    model_loc.freq = model_loc.freq(i_freq);
    if isfield(model, 'sigma')
        model_loc.sigma = model_loc.sigma(i_freq);
    end
    if isdistributed(Q) || iscodistributed(Q)
        Qloc = getLocalPart(Q);
        Qloc = Qloc(:,i_src,i_freq);
    else
        Qloc = Q(:,i_src);
    end
    Dloc = reshape(getLocalPart(D1),loc_size);
    if strcmp(dist_mode,'src')
        Dloc = permute(Dloc,[1 3 2]);
    end
    
    if params.wri > 0
        if length(params.lambda) > 1
            params.lambda = params.lambda(i_freq);
            params.i_freq = i_freq;
        end
    end
end

params_t = params{1};
if params_t.wri == 0
    clear params;
    params = params_t;
end

if obj_per_worker
    spmd,
        func = misfit_func_worker(Qloc,Dloc,model_loc,params,wri_mode);
    end
else
    if select_src_freq
        func = @(iS,iF) misfit_select_src_freq(iS,iF,i_src,i_freq,Qloc,Dloc,model_loc,params);
    else
        func = @(x) misfit_func(x,Qloc,Dloc,model_loc,params,wri_mode);
    end
end


end

function func = misfit_select_src_freq(iS,iF,i_src,i_freq,Qloc,Dloc,model_loc,params)

spmd,
    
end
func = @(x) misfit_func(x,Qloc,Dloc,model_loc,params);
end

function func = misfit_func_worker(Q,D,model,params,wri_mode)
func = @misfit_func1;
    function [f,g,h,w,f_aux] = misfit_func1(m)
        if isempty(Q) || isempty(D)
            f = 0; g = zeros(size(m)); h = 0; w = 0; f_aux = struct; f_aux.pde = 0; f_aux.dat = 0;
            return;
        end
        if wri_mode
            
            if params.FlagUandAlpha > 0
                [f,g,h,f_aux] = PDEfunc_wri_srcest2(PDEopts.OBJ,m,Q,[],D,model,params);
            else
                [f,g,h,f_aux] = PDEfunc_wri(PDEopts.OBJ,m,Q,[],D,model,params);
            end
        else
            [f,g] = PDEfunc(PDEopts.OBJ,m,Q,[],D,model,params);
            f_aux = struct;
            f_aux.pde = 0;
            f_aux.dat = f;
        end
        if wri_mode
            switch params.hessian
                case 'gn', h = opHGNpen(m,Q,D,model,params);
                case 'full', h = opHpen(m,Q,D,model,params);
            end
        else
            switch params.hessian
                case 'sparse', h = 0;
                case 'gn', h = opHGN(m,Q,model,params);
                case 'full', h = opH(m,Q,model,params);
            end
        end
        w = 0;
    end
end

function [f,g,h,w,f_aux] = misfit_func(m,Q,D,model,params,wri_mode)
spmd
    func = misfit_func_worker(Q,D,model,params,wri_mode);
    [f,g,h,w,f_aux] = func(m);
    f = pSPOT.utils.global_sum(f,1);
    g = pSPOT.utils.global_sum(g,1);
end
if isa(params,'Composite')
    p = params{1};
    hessian = p.hessian;
else
    hessian = params.hessian;
end
w = 0;
if wri_mode
    spmd,
        if strcmp(params.hessian,'sparse'), h = pSPOT.utils.global_sum(h,1); end
        f_dat = f_aux.dat; f_dat = pSPOT.utils.global_sum(f_dat,1);
        f_pde = f_aux.pde; f_pde = pSPOT.utils.global_sum(f_pde,1);
        if isfield(f_aux,'est_src')
            f_est_src = f_aux.est_src;
        else
            f_est_src = 0;
        end
    end
    C = params{1};
    if  C.FlagUandAlpha > 0
        
        %     g_est_src = f_aux.gest_src;
        g2 = [];
        
        for i = 1:length(f_aux)
            f_auxi = f_aux{i};
            g_est_src = f_auxi.gest_src;
            g2 = [g2; vec(g_est_src)];
            
        end
    end
    f_dat 	  = f_dat{1}; f_pde = f_pde{1};
    f_aux     = struct;
    f_aux.pde = f_pde;
    f_aux.dat = f_dat;
    
    f_aux.est_src = [];
    for i = 1: length(f_est_src)
        f_aux.est_src = [f_aux.est_src f_est_src{i}];
    end
    w = f_aux.est_src;
    
    switch hessian
        case 'sparse'
            h = h{1};
        case 'gn'
            h = oppHGNpen(m,Q,D,model,params);
        case 'full'
            h = oppHpen(m,Q,D,model,params);
        case 'hess-app'
            h = oppHpenApp(m,Q,D,model,params);
    end
else

    C = params;
    switch hessian
        case 'sparse'
            h = 0;
        case 'gn'
            h = oppHGN(m,Q,model,params);
        case 'full'
            h = oppH(m,Q,D,model,params);
    end
end
f  = f{1};
g  = g{1};
if  C.FlagUandAlpha > 0
    g  = [g;g2];
end


end
