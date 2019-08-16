function [ varargout ] = PDEfunc_dist( func, v, Q, input, Dobs,model,params )
%PDEFUNC_DIST - Parallel version of PDEfunc, distributed over (source x, source y, freq) pairs
%
% Curt Da Silva, 2015
% 
% Usage:
%   varargout = PDEfunc_dist( func, v, Q, input, Dobs, model, params);
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
%   params            - performance options struct
%      .srcfreqmask       - nsrc x nfreq binary matrix, an entry of 1 in the (i,j)th position indicates that this code should compute the desired quantity associated to the ith source and jth frequency (default: ones(nsrc,nfreq), all sources/freqs computed )
%
% Output:
%   if func=='obj',
%      {obj, gradient}
%
%   if func=='forw_model','jacob_forw'
%      forward, born-scattering wavefield, resp., of size [nrec x nsrc*nfreq], distributed over the last dimension
%
%   if func=='jacob_adj','hess_gn','hess'
%      migrated image, gauss-newton hessian, hessian image, resp., of size prod(model.n) x 1
%

% Parallelize over (source x, source y, frequency)
FORW_MODEL = PDEopts.FORW_MODEL; JACOB_FORW = PDEopts.JACOB_FORW; JACOB_ADJ = PDEopts.JACOB_ADJ;
HESS_GN = PDEopts.HESS_GN; HESS = PDEopts.HESS; OBJ = PDEopts.OBJ;
nfreq = length(model.freq);

if isfield(params,'wri'),wri_mode = params.wri; else wri_mode = false; end


if length(model.n)==2 || model.n(3)==1
    ndims = 2;
    nsrc = size(Q,2);
    nrec = length(model.zrec)*length(model.xrec);
    params.pdefunopts.numcompsrc = nsrc;
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
src_est_mode = params.pdefunopts.src_est_mode;
misfit_func = params.pdefunopts.misfit_func;

if ~isempty(params.pdefunopts.helm_dt) && norm(params.pdefunopts.helm_dt - model.d) > 0
    [~,~,subn] = fine2coarse(model.n,model.d,params.pdefunopts.helm_dt);            
    v = reshape(v,subn);
else
    v = reshape(v,model.n);
end
nlabs = parpool_size(); 

if ~isfield(params,'srcfreqmask')
    srcfreqmask = true(nsrc,nfreq);
else
    srcfreqmask = params.srcfreqmask;
    assert(size(srcfreqmask,1) == nsrc && size(srcfreqmask,2) == nfreq );
    assert(islogical(srcfreqmask));
end

nout = nargout;

% Serial version not recommended, you should do things in parallel for FWI
if nlabs == 0
    freqsxsy = cell(1,nfreq);
    for i=1:nfreq
        freqsxsy{i} = find(srcfreqmask(:,i));
    end
    if strcmp(src_est_mode,PDEopts.SRC_EST_NONE)~=0 && strcmp(func,OBJ)
        D = PDEfunc(FORW_MODEL,v,Q,[],[],model,params,freqsxsy);
        nfreq = length(find(any(srcfreqmask,1)));
        [Is,If] = ind2sub([nsrc,nfreq],find(vec(srcfreqmask)));
        If = If-find(any(srcfreqmask,1),1,'first')+1;
        w = est_src_weights(ones(nfreq,1),D,Dobs,[vec(Is),vec(If)],misfit_func);
        input = {w};
    end
    if nout==1
        output = PDEfunc(func,v,Q,input,Dobs,model,params,freqsxsy);
        varargout = {output};
    elseif nout==2
        [out1,out2] = PDEfunc(func,v,Q,input,Dobs,model,params,freqsxsy);
        varargout = {out1,out2};
    else
        [out1,out2,out3] = PDEfunc(func,v,Q,input,Dobs,model,params,freqsxsy);
        varargout = {out1,out2,out3};   
    end    
    return;
end

if ndims==2
    I = find(vec(srcfreqmask));
    minfreq = find(any(srcfreqmask,1),1,'first');
    npdes = length(I);
    spmd
        I = codistributed(I);  
        Dsize = [nsrc,nfreq];
        [isrc,ifreq] = ind2sub(Dsize,getLocalPart(I));
        
        freqs = unique(ifreq);
        freqsxsy = cell(1,nfreq);
        for i=1:length(freqs)
            freqsxsy{freqs(i)} = isrc(ifreq==freqs(i));
        end     
        ifreqs = find(~cellfun(@isempty,freqsxsy),1,'first'):find(~cellfun(@isempty,freqsxsy),1,'last');
        ifreqs = ifreqs-minfreq+1;
        if size(Q,3)>1 && (iscodistributed(Q)||isdistributed(Q))
            Q = getLocalPart(Q);
        end
    end    
else
    I = find(vec(srcfreqmask));
    npdes = length(I);
    % Setup distributed indices
    spmd
        I = codistributed(I);
        % Convert local 1D (src, freq) indices to cell array of per-frequency
        % source indices to compute, depending on which source dimension is
        % fixed
        switch fixed_src_dim
          case 1
            Dsize = [nsy,nsz,nfreq];
          case 2
            Dsize = [nsx,nsz,nfreq];
          case 3
            Dsize = [nsx,nsy,nfreq];
        end
        [isrc1,isrc2,ifreq] = ind2sub(Dsize,getLocalPart(I));    
        freqs = unique(ifreq);
        freqsxsy = cell(1,nfreq);
        for i=1:length(freqs)
            J = find(ifreq==freqs(i));
            freqsxsy{freqs(i)} = sub2ind(Dsize(1:2),vec(isrc1(J)),vec(isrc2(J)));
        end           
    end
end

% Source estimation
if strcmp(func,OBJ) && strcmp(src_est_mode,PDEopts.SRC_EST_NONE)==0
    if strcmp(src_est_mode,PDEopts.SRC_EST_RECOMPUTE)
        spmd
            D = PDEfunc(PDEopts.FORW_MODEL,v,Q,[],[],model,params,freqsxsy);
            D = codistributed.build(D,codistributor1d(2,codistributor1d.unsetPartition,[nrec,npdes]),'noCommunication');
        end        
        U = [];
    else
        spmd
            [D,U] = PDEfunc(PDEopts.FORW_MODEL,v,Q,[],[],model,params,freqsxsy);
            D = codistributed.build(D,codistributor1d(2,codistributor1d.unsetPartition,[nrec,npdes]),'noCommunication');
            U = codistributed.build(U,codistributor1d(2,codistributor1d.unsetPartition,[size(U,1),npdes]),'noCommunication');
        end        
    end
    nfreq = length(find(any(srcfreqmask,1)));
    [Is,If] = ind2sub([nsrc,nfreq],find(vec(srcfreqmask)));
    If = If-find(any(srcfreqmask,1),1,'first')+1;
    w = est_src_weights(ones(nfreq,1),D,Dobs,[vec(Is),vec(If)],misfit_func);
end

% Actual computational block
spmd
    if labindex~=1,params.disp_progress = false; params.pdefunopts.debug_mode = false; warning('off','all');  end 
    switch func
      case OBJ
        Dloc = getLocalPart(Dobs);
        if strcmp(src_est_mode,PDEopts.SRC_EST_RECOMPUTE)
            input = {w(ifreqs)};
        elseif strcmp(src_est_mode,PDEopts.SRC_EST_NORECOMPUTE)
            input = {w(ifreqs),getLocalPart(U)};
        else        
            input = {};
        end
        if nout == 1
            f = PDEfunc(func,v,Q,input,Dloc,model,params,freqsxsy);
        else
            [f,g,aux] = PDEfunc(func,v,Q,input,Dloc,model,params,freqsxsy);
            g = pSPOT.utils.global_sum(g,1);
            if wri_mode
                h = pSPOT.utils.global_sum(aux,1);                 
            else
                res = aux.res;
                res = codistributed.build(res,codistributor1d(1,codistributor1d.unsetPartition,[npdes,1]));
            end
        end            
        f = pSPOT.utils.global_sum(f,1); 
        dm = [];
      case JACOB_ADJ
        dm = PDEfunc(func,v,Q,getLocalPart(input),[],model,params,freqsxsy);
      case HESS
        dm = PDEfunc(func,v,Q,input,getLocalPart(Dobs),model,params,freqsxsy);
      case HESS_GN
        dm = PDEfunc(func,v,Q,input,[],model,params,freqsxsy);        
      case {FORW_MODEL,JACOB_FORW}
        D = PDEfunc(func, v, Q, input, [],model,params,freqsxsy);        
        dm = [];
        D = codistributed.build(D,codistributor1d(2,codistributor1d.unsetPartition,[nrec,npdes]),'noCommunication');
    end
    if ~isempty(dm), dm = pSPOT.utils.global_sum(dm,1); end
end



% Assign outputs
switch func
    case OBJ 
        if nout == 1
            varargout = {f{1}};
        else
            if wri_mode
                varargout = {f{1},g{1},h{1}};
            else
                aux = struct;
                res = gather(res);
                aux.res = res;
                if exist('w','var')
                    aux.src_weights = w;
                end
                varargout = {f{1},g{1},aux};
            end
        end
    case {JACOB_ADJ,HESS,HESS_GN}
        varargout = {dm{1}};
        
    case {FORW_MODEL,JACOB_FORW}
        varargout = {D};
end
            

end

