function [ varargout ] = PDEfunc_25D( func, v, Q, input, Dobs,model,params, freqsxsy )
% PDEfunc_25D - 2.5D PDEfunc
%
% Curt Da Silva, 2016
% 
% Usage:
%  varargout = PDEfunc_25D( func, v, Q, input, Dobs, model, params, freqsxsy );
%
% Input:
%   Same inputs as PDEfunc, except for 
% 
%   model.ysrc   - location of the source in the y-plane, single number
%   model.yrec   - location of the receivers in the y-plane, single number
% 
%   params
%      .nky      - number of 2D wavefields to compute (default: 100)
%      .use_knyq - if true, use the nyquist frequency (consistent gradients, but needs higher # of points) 
%                  if false, use a value slightly over the critical frequency (default)
%
     
OBJ = PDEopts.OBJ; FORW_MODEL = PDEopts.FORW_MODEL; JACOB_FORW = PDEopts.JACOB_FORW; JACOB_ADJ = PDEopts.JACOB_ADJ;

modes = {OBJ,FORW_MODEL,JACOB_FORW,JACOB_ADJ};

if ~any(ismember(func,modes)), error('Must choose a preset function.'); end
is_serial = @(x) ~isdistributed(x) && ~iscodistributed(x);
if ~isempty(Dobs), assert(is_serial(Dobs),'Need non-distributed data.'); end
if ~isempty(input), assert(is_serial(input),'Need non-distributed input.' ); end
assert(is_serial(Q), 'Need non-distributed source weights');
assert(isfield(model,'ysrc') && length(model.ysrc)==1,'length(model.ysrc) should be 1');
assert(isfield(model,'yrec') && length(model.yrec)==1,'length(model.yrec) should be 1');

lsopts = LinSolveOpts();
lsopts.solver = LinSolveOpts.SOLVE_LU;

params.scheme     = PDEopts.HELM2D_CHEN9P;
params.pml_max    = params.pdefunopts.helm_pml;
params.pml        = params.pdefunopts.helm_pml;
params.mat_free   = false;
params.solve_opts = lsopts;
if numel(v)==prod(model.n([1 3]))
    params.dt     = model.d([1 3]);
else
    params.dt     = params.pdefunopts.helm_dt;
end
if isfield(params,'use_knyq'),  use_knyq = params.use_knyq; else use_knyq = false; end
if isfield(params,'nky'), nky = params.nky; else nky = 100; end

nout = nargout;
assert(min(vec(v)) > 0,'Need a positive model parameter vector');

freq = model.freq;
nsrc = size(Q,2);
nrec = length(model.zrec)*length(model.xrec);
nfreq = length(freq);

% Deal with compatibility issues for old model structs
params = default_fwi_params2d(params);
model = fwi2d_model_compatibility(model);

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

if ~isempty(Dobs)
    Dobs = reshape(Dobs,nrec,npde_out);
end

srcw = fwi_wavelet(freq,model.t0,model.f0);

% Set outputs
switch func
  case OBJ
    f = 0;
    if nout >= 2
        g = zeros(numel(v),1);       
    end            
  case {FORW_MODEL,JACOB_FORW}   
    output = zeros(nrec,npde_out);        
  otherwise
    output = zeros(prod(model.n([1 3])),1);        
    if strcmp(func,JACOB_ADJ)
        input = reshape(input,nrec,npde_out);           
    end
end

switch model.unit
  case 'm/s', minv = min(vec(v));
  case 'km2/s2', minv = min(vec(1e3*(v.^(-1/2))));     
end

npdes = 1;

freq_idx = 0;
dy = model.d(2);
kynyq = pi/dy;
model_2d = model;
model_2d = rmfield(model_2d,'ysrc');
model_2d = rmfield(model_2d,'yrec');
% Permute so that the coordinates are (z,x), for compatibility with the 2D Helmholtz code
model_2d.n(2) = 1;
model_2d.n = model_2d.n([3 1 2]);
model_2d.d(2) = 1;
model_2d.d = model_2d.d([3 1 2]);
model_2d.o(2) = 0;
model_2d.o = model_2d.o([3 1 2]);

Hks = cell(nky,1);
u1 = cell(nky,1);
t = zeros(nky,1);
for k=1:nfreq
    if isempty(freqsxsy{k}); continue; end
    freq_idx = freq_idx+1;
    isrc = freqsxsy{k}; 
    
    freq = model.freq(k);
    
    % Critical frequency
    kc = 2*pi*freq/minv;
    
    % We're computing 1/pi \int_{a}^{b} u_{k_y}(z,x) cos(ky(y-y_src)) dky
    a = 0;
    
    % knyq will give you consistent gradients, but takes more points to 
    % get a nice looking wavefield
    if use_knyq       
        b = kynyq;
    else
        b = 1.2*kc;
    end
    c = 1/pi;
    
    % Gauss-Legendre quadrature
    [ky,w] = lgwt(nky,a,b);
    w = c*w;
        
    % Construct helmholtz matrices + quadrature weights
    for s=1:nky
        params.wn_offset = -(ky(s))^2;
    
        t(s) = w(s)*cos(ky(s)*(model.yrec-model.ysrc));        
        [Hk,comp_grid,dHk] = discrete_helmholtz(v,model_2d,freq,params);
        Hks{s} = Hk;
        
    end
    if strcmp(func,JACOB_FORW)
        dHdv = opMatrix(dHk)*opDiag_swp(comp_grid.phys2comp*vec(input));
    end
    
    % Source/receiver interpolation
    [zt,xt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);  
    Ps = opInterp('sinc',model.zsrc,zt,model.xsrc,xt);
    Pr = opInterp('sinc',model.zrec,zt,model.xrec,xt)';    
    
    current_src_idx = isrc;
    if length(size(Q))==3
        Qk_i = Ps*( srcw(k) * full(Q(:,current_src_idx,k)));
    else
        Qk_i = Ps*( srcw(k) * full(Q(:,current_src_idx)) );
    end
    
    % Scaling so that the wavefield amplitudes are the same for different grid spacings
    Qk_i = Qk_i * prod(model.d([1 3]))/prod(params.dt);       
    
    % Scaling to compensate for the y-direction
    Qk_i = Qk_i * dy;
    data_idx = npdes:npdes+length(current_src_idx)-1;
    
    u = zeros([prod(comp_grid.nt),length(current_src_idx)]);
    if strcmp(func,FORW_MODEL)|| strcmp(func,OBJ)
        % Solve for the field - weighted sum of 2D wavefields
        for s=1:nky            
            u1{s} = Hks{s}\Qk_i;
            u = u + u1{s}*t(s);
        end
        % Keep the wavefield signs consistent with the 3D Green's function
        u = conj(-u);
    end
    switch func
      case FORW_MODEL
        output(:,data_idx) = Pr*u;
      case OBJ
        r = Pr*u - Dobs(:,data_idx);
        f = f + 0.5*norm(r)^2;
        if nargout >=2 
            for s=1:nky
                v = Hks{s}'\(-Pr'*conj(-r));                
                output = output - t(s)*comp_grid.comp2phys*sum(real(conj(-u1{s}).* (dHk'*v)),2);
            end
        end
      case JACOB_FORW
        dout = 0;
        for s=1:nky
            u1{s} = Hks{s}\Qk_i;
            d = Hks{s}\(dHdv*(-u1{s}) );
            dout = dout + d*t(s);
        end
        dout = conj(-dout);
        output(:,data_idx) = Pr*dout;
      case JACOB_ADJ
        for s=1:nky
            u1{s} = Hks{s}\Qk_i;
            v = Hks{s}'\(-Pr'*conj(-input(:,data_idx)));                
            output = output - t(s)*comp_grid.comp2phys*sum(real(conj(-u1{s}) .* (dHk'*v)),2);
        end
    end
    npdes = npdes+length(current_src_idx);    
end

if strcmp(func,OBJ)
    if nargout==1
        varargout = {f};
    else
        varargout = {f,g};
    end
else
    if nout==1
        varargout = {output};
    else
        varargout = output; 
    end
end
    