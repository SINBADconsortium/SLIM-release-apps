function [ varargout ] = PDEfunc( func, v, Q, input, Dobs,model,params, freqsxsy )
%PDEFUNC3D - Abstract function for computing various quantities based on solutions of the 3D
%Helmholtz equation. Based on PDEfunc.m for the 2D case.
%
% Curt Da Silva, 2015
%
% Usage:
%   varargout = PDEfunc3D(func, v, Q, input, Dobs, model, params, freqsxsy);
%
% Input:
%   func               - string, one of
%                         - PDEopts.FORW_MODEL    - forward modelling operator
%                         - PDEopts.JACOB_FORW    - Jacobian of PDEopts.FORW_MODEL
%                         - PDEopts.JACOB_ADJ     - Jacobian adjoint of PDEopts.FORW_MODEL
%                         - PDEopts.HESS_GN       - Gauss-Newton hessian of PDEopts.FORW_MODEL
%                         - PDEopts.HESS          - Hessian of 'forw_model'
%                         - PDEopts.OBJ           - least-squares objective + gradient
%   v                 - vector with gridded velocity parameter (units
%                       defined by model.unit)
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source
%   input             - input model perturbation, used by 'jacob_forw','hess_gn','hess'
%                       or data residual used by 'jacob_adj'
%   Dobs              - observed data, used by 'hess','obj'
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.freq        - frequencies (in Hz)
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc,ysrc} - vectors describing source array
%   model.{zrec,xrec,yrec} - vectors describing receiver array.
%
%   params            -
%     .lsopts         - a LinSolveOpts object, specifying the parameters used for
%                       solving the helmholtz equations
%     .pdefunopts     - a PDEopts object, specifying the parameters used in
%                       this function
%
%   freqsxsy          - nfreq x 1 cell array, each entry is a nsrc_comp x 1
%                       array of (1D) indices or a nsrc_comp x 2 array of (2D) indices
%                       that indicates what sources should be computed for this frequency
%                       (default: )
%
% Output:
%   if func=='obj',
%      {obj, gradient}
%
%   if func=='forw_model','jacob_forw'
%      forward, born-scattering wavefield, resp., of size [nrec x nsrc*nfreq]
%
%   if func=='jacob_adj','hess_gn','hess'
%      migrated image, gauss-newton hessian, hessian image, resp., of size prod(model.n) x 1
%
%

FORW_MODEL = PDEopts.FORW_MODEL; JACOB_FORW = PDEopts.JACOB_FORW; JACOB_ADJ = PDEopts.JACOB_ADJ;
HESS_GN = PDEopts.HESS_GN; HESS = PDEopts.HESS; OBJ = PDEopts.OBJ; FIELD = PDEopts.FIELD;

modes = {OBJ,FORW_MODEL,JACOB_FORW,JACOB_ADJ,HESS_GN,HESS};

nout = nargout;

% Sanity checks
if ~any(ismember(func,modes)), error('Must choose a preset function.'); end
if ~isempty(input), assert(~isdistributed(input) && ~iscodistributed(input),'Need non-distributed input.' ); end
assert(~isdistributed(Q) && ~iscodistributed(Q), 'Need non-distributed source weights');

% Unpack options, reformulate them for further propagation
if ~isfield(params,'pdefunopts') || ~isa(params.pdefunopts,'PDEopts')    
    error('Must specify pdefuncopts of type PDEopts');
else
    % parameters that are used by this function
    zero_boundary       = params.pdefunopts.zero_boundary;
    window_source_grad  = params.pdefunopts.window_source_grad;
    src_est_mode        = params.pdefunopts.src_est_mode;
    rec_interp          = params.pdefunopts.rec_interp;
    src_interp          = params.pdefunopts.src_interp;
    debug_mode          = params.pdefunopts.debug_mode;
    misfit              = params.pdefunopts.misfit_func;
    if isfield(params,'hessian'), hessian_mode = params.hessian; else hessian_mode = PDEopts.HESS_NONE; end
    if isfield(params,'wri'), wri_mode = params.wri; else wri_mode = false; end
    if wri_mode, lambdas = params.lambda; end
    
    % parameters that get propagated to lower level functions
    params.scheme       = params.pdefunopts.helm_scheme;
    if numel(v)==prod(model.n)
        params.dt       = model.d;
    else
        params.dt       = params.pdefunopts.helm_dt;
    end
    params.free_surface = params.pdefunopts.helm_free_surface;
    params.pml_max      = params.pdefunopts.helm_pml_max;    
    params.pml          = params.pdefunopts.helm_pml;
    params.mat_free     = params.pdefunopts.helm_mat_free;
    params.n_threads    = params.pdefunopts.helm_nthreads;
    params.disp_output  = debug_mode;    
end


if ~isfield(params,'lsopts') || ~isa(params.lsopts,'LinSolveOpts')
    error('Must specify lsopts of type LinSolveOpts');
else
    params.solve_opts   = params.lsopts;   
end

if isfield(params,'disp_progress'),disp_progress = params.disp_progress; else disp_progress = false; end

assert(min(vec(v)) > 0,'Need a positive model parameter vector');

% Extract geometry information + what sources/frequencies are computed
freq = model.freq;
nfreq = length(freq);

if ~isempty(Dobs)
    if isa(Dobs,'double')
        if isfield(model,'srcgrid') %irregular source grid
            Dobs = FWIFreqData(Dobs,model.srcgrid,model.freq);
        else %regular source grid            
            Dobs = FWIFreqData(Dobs,model);
        end
    else
        assert(isa(Dobs,'FWIFreqData'));
    end
end

if length(model.n)==2 ||  model.n(3)==1
    ndims = 2; 
    params = default_fwi_params2d(params);
    model = fwi2d_model_compatibility(model);
    if isempty(Dobs)        
        if ~isfield(model,'srcgrid')
            sg = RegularGrid(model.zsrc,model.xsrc);
            rg = RegularGrid(model.zrec,model.xrec);
            srcgrid = DataGrid(sg,rg);
        else
            srcgrid = model.srcgrid;
        end
    else
        srcgrid = Dobs.getSourceGrid();
    end
else
    ndims = 3;
    if ~isempty(Dobs) 
        srcgrid = Dobs.getSourceGrid();
    elseif isfield(model,'srcgrid') && isa(model,'DataGrid')
        srcgrid = model.srcgrid;
    else
        sg = RegularGrid(model.xsrc,model.ysrc,model.zsrc);
        rg = RegularGrid(model.xrec,model.yrec,model.zrec);
        srcgrid = DataGrid(sg,rg);
    end
    
end
nsrc = srcgrid.numsrcs();    
nrec = srcgrid.numrecs();
if ndims==2, numcompsrc = nsrc; else numcompsrc = 1; end

% If unspecified for 2D, assume the model vector is s2/km2, for
% legacy examples
if ndims==2 && ~isfield(model,'unit'), model.unit = 's2/km2'; end

% Parse freqsxsy matrix, which is a logical nsrc x nfreq matrix denoting which (src,freq) pairs to compute, in to a nfreq x 1 cell array containing those indices
assert(exist('freqsxsy','var')==0 || isempty(freqsxsy) || (islogical(freqsxsy) && numel(freqsxsy )==nsrc*nfreq),'freqsxsy should be logical and of size nsrc x nfreq');

if exist('freqsxsy','var')==0 || isempty(freqsxsy)
    freqsxsy = cell(1,nfreq);
    for i=1:length(freqsxsy)
        freqsxsy{i} = vec(1:nsrc);
    end
    npde_out = nsrc*nfreq;
else
    freqsxsy = reshape(freqsxsy,nsrc,nfreq);
    if ndims==3 && length(find(any(freqsxsy,1))) > 1
        error('Only a single frequency at a time supported for 3D FWI');
    end
    I = find(vec(freqsxsy));    
    npde_out = length(I);
    [iS,iF] = ind2sub([nsrc,nfreq],I);
    freqsxsy = cell(1,nfreq);
    for i=1:length(freqsxsy)
        J = find(iF == i);
        if ~isempty(J)
            freqsxsy{i} = iS(J);
        end
    end           
end
% If we don't have any work to do, return the right outputs (zeros, of the right size)
if npde_out==0
    z = zeros(prod(model.n),1);
    switch func
      case OBJ
        if nout >= 3
            varargout = {0,z,z};            
        elseif nout >= 2
            varargout = {0,z};
        else
            varargout = {0};
        end
      case {FORW_MODEL,JACOB_FORW}
        varargout = {zeros(nrec,0)};      
      otherwise
        varargout = {z};        
    end
    return;
end


if model.f0==0
    error('Must have a positive peak frequency');
end

% Define wavelet
w = fwi_wavelet(freq,model.t0,model.f0);
U_computed = []; 
src_weights = ones(nfreq,1);
if ~strcmp(src_est_mode,PDEopts.SRC_EST_NONE) && ~isempty(input)
    weights = input{1};
    src_weights = weights;
    if length(input)>1
        U_computed = input{2};
    end
end

% Set outputs
switch func
  case OBJ
    f = 0;
    if nout >= 2
        g = zeros(numel(v),1);
        if wri_mode
            h = sparse([],[],[],prod(model.n),prod(model.n)); 
        else
            h = zeros(numel(v),1);
        end
    end            
  case {FORW_MODEL,JACOB_FORW}   
    output = zeros(nrec,npde_out);        
  otherwise
    output = zeros(prod(model.n),1);        
    if strcmp(func,JACOB_ADJ)
        input = reshape(input,nrec,npde_out);           
    end
end


% Main computational loop
npdes = 1;

% progress indicator variables
disp_out_freq = 0.25:0.25:1; disp_counter = 1; s = tic;
freq_idx = 0;

for k=1:nfreq
    if isempty(freqsxsy{k}), continue; end
    freq_idx = freq_idx+1;
    isrc = freqsxsy{k}; 
            
    % Chunk up the sources in isrc in to blocks of size numcompsrc
    src_batches = index_block(isrc,numcompsrc);
    
    if wri_mode
        if length(lambdas)>1
            lambda = lambdas(k);
        else
            lambda = lambdas;
        end
    end 
    
    % We have several grids in play here
    % 1) Model grid : (possible coarsening of) the finest grid on which the model lives
    % 2) Pml grid: model grid + pml padding
    
    % 3) Source grid : grid on which the sources reside
    % 4) Receiver grid : grid on which the receivers reside
    
    
    % Get helmholtz matrix + computational grid parameters
    % The grid spacing parameter is specified by calling function
    % H : pml grid -> pml grid
    if strcmp(func,HESS) || wri_mode
        [H,comp_grid,T,DT_adj] = discrete_helmholtz(v,model,freq(k),params);        
    else
        [H,comp_grid,T] = discrete_helmholtz(v,model,freq(k),params);
        ddH = [];
    end    
    
    if npdes==1 && nout>=2 && strcmp(func,FORW_MODEL)
        fields = zeros(prod(comp_grid.nt),npde_out);    
    else
        fields = zeros(prod(comp_grid.nt),0);
    end
    
    % operators from pml grid -> model grid + vice versa
    to_comp = comp_grid.phys2comp;
    to_phys = comp_grid.comp2phys;
    
    % WRI specific operators
    if wri_mode
        [Ps,Pr] = src_rec_interp(ndims,model,comp_grid,src_interp,rec_interp);
        H = H.params.coef;
        Pe = comp_grid.comp2phys';
        
        % Indicator function for the computational grid in the extended domain
        idx_m = logical( Pe*(ones(size(Pe,2),1)) );
        
        e = zeros(size(Pr,1),1); e(1) = 1;
        I = []; J = []; V = [];
        z = zeros(size(Pr,1),1);
        for i=1:size(Pr,1)
            if i > 1, e(i-1) = 0; end, e(i) = 1;
            z = Pr'*e;
            J1 = find(z);
            I = [I;i*ones(length(J1),1)]; J = [J;vec(J1)]; V = [V;vec(z(J1))];    
        end
        Prsp = sparse(I,J,V,size(Pr,1),size(Pr,2));
    
        if ~isempty(input), dm1 = comp_grid.phys2comp*input; 
        else dm1 = []; end
        A = [lambda*H; Prsp];
    end            
    
    % Loop over all source batches
    for j=1:length(src_batches)        
        current_src_idx = src_batches{j};                        
                
        % Inject source weights in to the computational grid
        if length(size(Q))==3
            active_srcs = find(any(Q(:,current_src_idx,k),2));
            [Ps,Pr] = srcgrid.interp_operators(active_srcs,comp_grid,src_interp);
            Qk_i = Ps*( w(k) * full(Q(:,current_src_idx,k)));
        else
            active_srcs = find(any(Q(:,current_src_idx),2));
            [Ps,Pr] = srcgrid.interp_operators(active_srcs,comp_grid,src_interp);
            Qk_i = Ps*( w(k)* full(Q(:,current_src_idx)) );
        end
        
        % Scaling so that the wavefield amplitudes are the same for different grid spacings
        Qk_i = Qk_i * prod(model.d(1:ndims))/prod(params.dt(1:ndims));       
                
        data_idx = npdes:npdes+length(current_src_idx)-1;
        if isempty(U_computed) 
            % Primary wavefield - prod(comp_grid.nt) x ncompsrc matrix
            if wri_mode
                U = [];
            else
                U = H \ Qk_i;
            end
        else
            U = U_computed(:,data_idx);
        end
                
        if debug_mode
            assert(norm(H*U-Qk_i,'fro') <= 10* params.lsopts.tol * norm(Qk_i,'fro')*numcompsrc);
        end                       
        
        if wri_mode                        
            if nargout==1
                out = PDEfunc_compute_wri(func,A,getData(Dobs,data_idx),U,Qk_i,src_weights(freq_idx),Pr,H,dH_forw,T_adj,lambda,dH,ddH,dm1,to_phys,idx_m);
            elseif nargout==2
                [f1,g1] = PDEfunc_compute_wri(func,A,getData(Dobs,data_idx),U,Qk_i,src_weights(freq_idx),Pr,H,lambda,dH,ddH,dm1,to_phys,idx_m);
            elseif nargout==3
                [f1,g1,h1] = PDEfunc_compute_wri(func,A,getData(Dobs,data_idx),U,Qk_i,src_weights(freq_idx),Pr,H,lambda,dH,ddH,dm1,to_phys,idx_m);
            end
            switch func
              case FIELD
                output(:,current_src_idx,k)=out;
              case FORW_MODEL
                output(:,data_idx) = out;
              case OBJ
                if nargout==1
                    f = f + out;
                elseif nargout>=2
                    f = f + f1;
                    g = g + g1;                   
                    if nargout==3
                        h = h + h1;
                    end
                end
              case {HESS_GN,HESS}
                output = output + out;
            end
        else
            % Apply source weights
            U = U*src_weights(freq_idx);
            sum_srcs = @(x) to_phys*sum(real(x),2);
            switch func
              case OBJ            
                [phi,dphi] = misfit(Pr*U,getData(Dobs,data_idx),current_src_idx,freq_idx);
                
                f = f + phi;
                if nargout >= 2
                    V = H' \ ( -Pr'* dphi);
                    g = g + sum_srcs(T(U)'*V);
                    if debug_mode
                        assert(norm(H'*V-(-Pr'*dphi),'fro') <= 10* params.lsopts.tol * norm(-Pr'*dphi,'fro')*numcompsrc);
                    end                    
                    
                    if nargout >= 3      
                        switch hessian_mode
                          case PDEopts.HESS_GN_DIAG
                            % True Gauss-Newton diagonal
                            W = H'\(Pr'*eye(size(Pr,1)));                            
                            h = h + to_phys* sum(repmat(vec(sum( abs(W'*dH) .^2,1)),1,size(U,2)) .* (abs(U).^2),2);
                            
                          case PDEopts.HESS_DIAG_SHIN01
                            % Pseudo hessian computation of Shin 2001
                            h = h + to_phys* sum(repmat(vec(sum((Pr*dH).^2,1)),1,size(U,2)).* (abs(U).^2),2);
                            
                          case PDEopts.HESS_DIAG_ENCODE
                            % number of sim receivers
                            m = 3;
                            S = randn(size(Pr,1),m)/(sqrt(m));
                            W = H'\( Pr'*S );                            
                            h = h + to_phys* sum(repmat(vec(sum( abs(W'*dH) .^2,1)),1,size(U,2)) .* (abs(U).^2),2);
                          otherwise
                            h = 0;
                        end
                            
                    end
                end            
              case FORW_MODEL
                output(:,data_idx) = Pr*U;
                if size(fields,2)>0, fields(:,data_idx) = U; end
                
              case JACOB_FORW
                dm = to_comp*vec(input);
                dU = H\(-T(U)*dm);
                output(:,data_idx) = Pr*dU;
                
              case JACOB_ADJ
                V = H'\( -Pr'* input(:,data_idx) );                
                output = output + sum_srcs( T(U)'*V );
                if debug_mode
                    assert(norm(H'*V-(-Pr'*input(:,data_idx))) <  params.lsopts.tol * norm(input(:,data_idx))*numcompsrc);
                end
                
              case HESS_GN
                dm = to_comp*vec(input);
                dU = H\(-T(U)*dm);
                dU = H'\(-Pr'*Pr*dU);
                output = output + sum_srcs(T(U)'*dU);
                
              case HESS
                dm = to_comp*vec(input);
                [~,dphi,d2phi] = misfit(Pr*U,getData(Dobs,data_idx),current_src_idx,freq_idx);
                dU = H\(-T(U)*dm);                            
                V = H'\(-Pr'*dphi);
                dV = H'\(-T(V)*dm - Pr'* reshape(d2phi*vec(Pr*dU),nrec,size(U,2)));
                output = output + sum_srcs(DT_adj(U,dm,dU)*V + T(U)'*dV);
            end
        end
        npdes = npdes+length(current_src_idx);
        
        if disp_progress, 
            if disp_counter < length(disp_out_freq) && (npdes-1)/npde_out >= disp_out_freq(disp_counter)
                s = toc(s);
                work_left = (1-(npdes-1)/npde_out)/max(diff(disp_out_freq));
                secs_left = work_left*s;
                disp([datestr(now) ' : ' num2str((npdes-1)/npde_out*100) '% finished, ETA ' datestr(now + seconds(secs_left)) ]);
                disp_counter = disp_counter+1;
                s = tic;
            end
        end
    end    
end
if disp_progress
    disp([datestr(now) ' : 100% finished ']);
end

if ~isempty(fields) && numel(fields) > 0, output = {output,fields}; end

% Set outputs and window, if specified by user
if strcmp(func,OBJ)
    if nargout == 1
        varargout = {f};
    else
        if zero_boundary,
            g = reshape(g,model.n);
            g([1 end],:,:) = 0; g(:,[1 end],:) = 0; g(:,:,[1 end]) = 0;
            g = vec(g);
        end
        if window_source_grad
            M = src_mask(model,2*model.d(1));
            g = g.*vec(M);
        end
        if wri_mode
            varargout = {f,g,h};
        else
            varargout = {f,g,h};
        end
    end
else    
    if nout==1
        varargout = {output};
    else
        varargout = output;
    end
end


end
