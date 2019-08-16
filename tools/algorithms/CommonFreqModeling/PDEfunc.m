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

SRC_INTERP_LIN = PDEopts.SRC_INTERP_LIN; SRC_INTERP_SINC = PDEopts.SRC_INTERP_SINC;
REC_INTERP_LIN = PDEopts.REC_INTERP_LIN; REC_INTERP_SINC = PDEopts.REC_INTERP_SINC;

nout = nargout;

% Sanity checks
if ~any(ismember(func,modes)), error('Must choose a preset function.'); end
is_serial = @(x) ~isdistributed(x) && ~iscodistributed(x);
if ~isempty(Dobs), assert(is_serial(Dobs),'Need non-distributed data.'); end
if ~isempty(input), assert(is_serial(input),'Need non-distributed input.' ); end
assert(is_serial(Q), 'Need non-distributed source weights');

% Unpack options, reformulate them for further propagation
if ~isfield(params,'pdefunopts') || ~isa(params.pdefunopts,'PDEopts')    
    error('Must specify pdefuncopts of type PDEopts');
else
    % parameters that are used by this function
    numcompsrc          = params.pdefunopts.numcompsrc;
    zero_boundary       = params.pdefunopts.zero_boundary;
    window_source_grad  = params.pdefunopts.window_source_grad;
    src_est_mode        = params.pdefunopts.src_est_mode;
    rec_interp          = params.pdefunopts.rec_interp;
    src_interp          = params.pdefunopts.src_interp;
    debug_mode          = params.pdefunopts.debug_mode;
    misfit              = params.pdefunopts.misfit_func;
    offset_mask         = params.pdefunopts.offset_mask;
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

if length(model.n)==2 ||  model.n(3)==1
    ndims = 2;
    nsrc = size(Q,2);
    nrec = length(model.zrec)*length(model.xrec);
    numcompsrc = nsrc;
    params = default_fwi_params2d(params);
    model = fwi2d_model_compatibility(model);

else
    ndims = 3;
    nsx = length(model.xsrc); nsy = length(model.ysrc); 

    nsz = length(model.zsrc);
    if nsx == 1, fixed_src_dim = 1; elseif nsy==1, fixed_src_dim = 2; elseif nsz==1, fixed_src_dim = 3; end
    switch fixed_src_dim
      case 1
        src_dims = [nsy,nsz];
      case 2
        src_dims = [nsx,nsz];
      case 3
        src_dims = [nsx,nsy];
    end
    nsrc = nsx*nsy*nsz;
    nrec = length(model.yrec)*length(model.xrec)*length(model.zrec);
end
% If unspecified for 2D, assume the model vector is s2/km2, for
% legacy examples
if ndims==2 && ~isfield(model,'unit'), model.unit = 's2/km2'; end
nfreq = length(freq);
if isempty(offset_mask), offset_mask = ones(nrec,nsrc); end

if exist('freqsxsy','var')==0
    freqsxsy = cell(1,nfreq);
    for i=1:length(freqsxsy)
        freqsxsy{i} = vec(1:nsrc);
    end
    npde_out = nsrc*nfreq;
else
    assert(length(freqsxsy)==nfreq);
    npde_out = 0;
    for i=1:length(freqsxsy)
        if ~isempty(freqsxsy{i})
            if min(size(freqsxsy))> 1
                % 2D indices, convert to 1D indices                
                I = freqsxsy{i};
                I = ind2sub(src_dims,I(:,1),I(:,2));
                freqsxsy{i} = I;                           
            end
            npde_out = npde_out + numel(freqsxsy{i});
        end
    end
end
% If we don't have any work to do, return the right outputs (zeros, of the right size)
if npde_out==0
    z = zeros(prod(model.n),1);
    switch func
      case OBJ
        if nout >= 3
            s = struct;
            s.res = zeros(0,1);
            varargout = {0,z,s};            
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

if ~isempty(Dobs)
    Dobs = reshape(Dobs,nrec,npde_out);
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

% Auxiliary outputs
if strcmp(func,OBJ)
    aux = struct;
    aux.res = zeros(npde_out,1);
end


% Main computational loop
npdes = 1;

% progress indicator variables
disp_out_freq = 0.25:0.25:1; disp_counter = 1; T = tic;
freq_idx = 0;
fields = [];
for k=1:nfreq
    if isempty(freqsxsy{k}), continue; end
    freq_idx = freq_idx+1;
    isrc = freqsxsy{k}; nsrcloc = length(isrc);
    
    % Helmholtz solves are performed in batches of size determined by numcompsrc
    % Make sure that the Matlab worker has sufficient memory/is efficient enough to process
    % numcompsrc sources in parallel
    if nsrcloc <= numcompsrc,
        src_batches = cell(1,1); src_batches{1} = isrc;
    else
        kup = floor(nsrcloc/numcompsrc);
        src_batches = cell(kup,1);
        src_total = 1;
        for j=1:kup
            src_batches{j} = isrc(src_total:src_total+numcompsrc-1);
            src_total = src_total + numcompsrc;
        end
        if kup*numcompsrc < nsrcloc,
            if nsrcloc-kup*numcompsrc < 0.5*numcompsrc
                %Too few extra sources, distribute them to the last round
                src_batches{end} = [src_batches{end};isrc(src_total:nsrcloc)];
            else
                %Number of left over sources is efficient enough to have its own computational round
                src_batches{end+1} = src_total:nsrcloc;
            end
        end
    end
    
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
    % Hk : pml grid -> pml grid
    if strcmp(func,HESS) || wri_mode
        [Hk,comp_grid,dH,ddH] = discrete_helmholtz(v,model,freq(k),params);        
    else
        [Hk,comp_grid,dH] = discrete_helmholtz(v,model,freq(k),params);
    end
    
    
    if npdes==1 && nout>=2 && strcmp(func,FORW_MODEL)
        fields = zeros(prod(comp_grid.nt),npde_out);    
    end
    
    % operators from pml grid -> model grid + vice versa
    to_comp = comp_grid.phys2comp;
    to_phys = comp_grid.comp2phys;
    if strcmp(func,JACOB_FORW) || strcmp(func,HESS_GN) || strcmp(func,HESS)
        if ndims==2
            dHdv = opMatrix(dH) * opDiag_swp(to_comp*vec(input));
            if strcmp(func,HESS), ddHdv = opMatrix(ddH)*opDiag_swp(to_comp*vec(input)); end
        else
            dHdv = dH*opDiag(to_comp*vec(input));
            if strcmp(func,HESS), ddHdv = ddH*opDiag(to_comp*vec(input)); end
        end
    end
    
    % Source / receiver interpolation
    if ndims==2
        % 2D model is organized as (z,x) for legacy/visualization reasons
        [zt,xt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);           
        if strcmp(src_interp,SRC_INTERP_LIN)
            Ps = opInterp('cubic',zt,model.zsrc,xt,model.xsrc)';            
        elseif strcmp(src_interp,SRC_INTERP_SINC)
            Ps = opInterp('sinc',model.zsrc,zt,model.xsrc,xt);
        end
        if strcmp(rec_interp,REC_INTERP_LIN)
            Pr = opInterp('cubic',zt,model.zrec,xt,model.xrec);
        elseif strcmp(rec_interp,REC_INTERP_SINC)
            Pr = opInterp('sinc',model.zrec,zt,model.xrec,xt)'; 
        end
    else
        % 3D model is (x,y,z)
        [xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
        if strcmp(src_interp,SRC_INTERP_LIN)
            Ps = opInterp('cubic',xt,model.xsrc,yt,model.ysrc,zt,model.zsrc)';     
        elseif strcmp(src_interp,SRC_INTERP_SINC)
            Ps = opInterp('sinc',model.xsrc,xt,model.ysrc,yt,model.zsrc,zt);     
        end
        if strcmp(rec_interp,REC_INTERP_LIN)
            Pr = opInterp('cubic',zt,model.zrec,xt,model.xrec);            
        elseif strcmp(rec_interp,REC_INTERP_SINC)
            Pr = opInterp('sinc',model.xrec,xt,model.yrec,yt,model.zrec,zt)';
        end        
    end           
    
    % WRI specific operators
    if wri_mode
        Hk = Hk.params.coef;
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
    
        if ~isempty(input), dm1 = comp_grid.phys2comp*input; end
        A = [lambda*Hk; Prsp];
    end
    
    
    % Loop over all source batches
    for j=1:length(src_batches)        
        current_src_idx = src_batches{j};
        W = offset_mask(:,current_src_idx);
        
        % Inject source weights in to the computational grid
        if length(size(Q))==3
            Qk_i = Ps*( w(k) * full(Q(:,current_src_idx,k)));
        else
            Qk_i = Ps*( w(k)* full(Q(:,current_src_idx)) );
        end
        
        % Scaling so that the wavefield amplitudes are the same for different grid spacings
        Qk_i = Qk_i * prod(model.d(1:ndims))/prod(params.dt(1:ndims));       
        
        
        data_idx = npdes:npdes+length(current_src_idx)-1;
        if isempty(U_computed)
            % Primary wavefield - prod(comp_grid.nt) x ncompsrc matrix
            if wri_mode
                S = [lambda*Qk_i;Dobs(:,data_idx)];
                G = A'*A;
                Uk = A\S;
                Dp = Pr*Uk;
                pde_res = Hk*Uk-Qk_i;
            else
                Uk = Hk \ Qk_i;
            end
        else
            Uk = U_computed(:,data_idx);
        end
                
        if debug_mode
            assert(norm(Hk*Uk-Qk_i,'fro') <= 10* params.lsopts.tol * norm(Qk_i,'fro')*numcompsrc);
        end       
        
        % Apply source weights
        Uk = Uk*src_weights(freq_idx);
        
        if wri_mode
            switch func
              case FIELD
                output(:,current_src_idx,k)=Uk;
              case FORW_MODEL
                output(:,data_idx) = Dp;
              case OBJ
                fdata = 0.5 * norm(W.*(Dp-Dobs(:,data_idx)));
                fpde = 0.5* lambda^2*norm(pde_res,'fro')^2;
                f = f + fdata + fpde;
                if nargout >= 2
                    V1 = conj(Uk); V2 = dH'*pde_res;
                    g = g + lambda^2 * to_phys * (sum(real((V1.*V2)),2));
                    ddH1   = ddH(idx_m,idx_m);
                    dH1    = dH(idx_m,idx_m);
                    for i=1:size(Uk,2)          
                        U_diag  = Uk(:,i); U_diag = spdiags(U_diag(idx_m),0,prod(model.n),prod(model.n));
                        res_diag = pde_res(:,i); res_diag = spdiags(res_diag(idx_m),0,prod(model.n),prod(model.n));
                        O   = ddH1'*res_diag;
                        P   = dH1*U_diag;                              
                        h   = h + lambda^2 * real( U_diag'*O +  P'*P);
                    end
                end
              case HESS_GN
                dHdm   = dH*opDiag_swp(dm1);
                dUk     = lambda^2 * ( G \ (-dHdm'*(pde_res) - Hk' * (dHdm*Uk)));
                Z1     = Pr*dUk; Z2 = lambda * (dHdm * Uk + Hk * dUk);
                Z1     = Pr' * Z1 + lambda * Hk' * Z2;
                Z1     = lambda^2 * ( G \ Z1 );
                R      = conj( -dH' * (pde_res)) .* Z1;
                R      = sum(R - conj(Uk)  .* (dH' * Hk * Z1),2);
                Z2     = sum( conj( Uk ) .* (dH' * (lambda*Z2)) , 2);
                output = output + to_phys*(real(R + Z2));
              case HESS
                dHdm   = dH*opDiag_swp(dm1);
                ddHtdm  = ddH'*opDiag_swp(dm1);
                dUk     = lambda^2 * (G \ (-dHdm'*(pde_res) - Hk' * (dHdm*Uk)));
                % d/dm^2 phi(m,u)[dm] 
                output = output + lambda^2*to_phys*sum(real( conj(Uk) .* (dH'*(dHdm*Uk)) ),2);
                % d/dm d/du phi(m,u) [du]
                output = output + lambda^2*to_phys*sum(real( conj(Uk) .* (ddHtdm*pde_res)),2);
                output = output + lambda^2*to_phys*sum(real( conj(dUk) .* (dH' *pde_res)  ),2);
                output = output + lambda^2*to_phys*sum(real( conj(Uk) .* ( dH' *(Hk*dUk) ) ),2);   
            end
        else
            switch func
              case OBJ            
                [phi,dphi] = misfit(W.*(Pr*Uk),W.*squeeze(Dobs(:,data_idx)),current_src_idx,freq_idx);
                
                f = f + phi;
                if nargout >= 2
                    Vk = Hk' \ ( -Pr'* dphi);
                    g = g + to_phys*sum(real(conj(Uk) .* (dH'*Vk)),2);
                    if debug_mode
                        assert(norm(Hk'*Vk-(-Pr'*r),'fro') <= 10* params.lsopts.tol * norm(-Pr'*r,'fro')*numcompsrc);
                    end
                end            
              case FORW_MODEL
                output(:,data_idx) = Pr*Uk;
                if ~isempty(fields), fields(:,data_idx) = Uk; end
              case JACOB_FORW
                dUk = Hk\(dHdv*(-Uk));
                output(:,data_idx) = Pr*dUk;
                
              case JACOB_ADJ
                Vk = Hk'\( -Pr'* input(:,data_idx) );
                output = output + to_phys*sum(real(conj(Uk) .* (dH'*Vk)),2);
                if debug_mode
                    assert(norm(Hk'*Vk-(-Pr'*input(:,data_idx))) <  params.lsopts.tol * norm(input(:,data_idx))*numcompsrc);
                end
                
              case HESS_GN
                dUk = Hk\(-dHdv*Uk);
                dUk = Hk'\(-Pr'*Pr*dUk);
                output = output + to_phys*sum(real(conj(Uk) .* (dH'*dUk)),2);
                
              case HESS
                [~,dphi,d2phi] = misfit(Pr*Uk,squeeze(Dobs(:,data_idx)),current_src_idx,freq_idx);            
                dUk = Hk\(-dHdv*(Uk));            
                V1  = Uk;                 dV1 = dUk;
                V2  = dH';                dV2 = ddHdv';
                V3  = Hk'\( -Pr'* dphi);  dV3 = Hk'\ ( -dHdv'* V3  - Pr'*d2phi*Pr*dUk );
                
                output = output + to_phys*sum(real(conj(dV1) .* (V2 * V3) + ...
                                                   conj(V1) .* (dV2*V3 + V2*dV3)),2);
                
            end
        end
        npdes = npdes+length(current_src_idx);
        
        if disp_progress, 
            if disp_counter < length(disp_out_freq) && (npdes-1)/npde_out >= disp_out_freq(disp_counter)
                T = toc(T);
                work_left = (1-(npdes-1)/npde_out)/max(diff(disp_out_freq));
                secs_left = work_left*T;
                disp([datestr(now) ' : ' num2str((npdes-1)/npde_out*100) '% finished, ETA ' datestr(now + seconds(secs_left)) ]);
                disp_counter = disp_counter+1;
                T = tic;
            end
        end
    end
end

if ~isempty(fields), output = {output,fields}; end

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
            varargout = {f,g,aux};
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
