function [ varargout ] = PDEfunc_dist( func, v, Q, input, Dobs,model,params,srcfreqmask )
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
if isfield(params,'disp_host'),disp_host = params.disp_host; else disp_host = false; end
if isfield(params,'wri'),wri_mode = params.wri; else wri_mode = false; end

nsrc = size(Q,2);
if length(model.n)==2 || model.n(3)==1
    ndims = 2;
    nrec = length(model.zrec)*length(model.xrec);
    params.pdefunopts.numcompsrc = nsrc;
    params = default_fwi_params2d(params);
    model = fwi2d_model_compatibility(model);
else
    ndims = 3;
    nrec = length(model.yrec)*length(model.xrec)*length(model.zrec);
end
src_est_mode = params.pdefunopts.src_est_mode;
misfit_func = params.pdefunopts.misfit_func;

% If the user specifies a dt for the velocity (vector) which is different from model.d, reshape it 
if ~isempty(params.pdefunopts.helm_dt) && norm(params.pdefunopts.helm_dt - model.d) > 0
    [~,~,subn] = fine2coarse(model.n,model.d,params.pdefunopts.helm_dt);            
    v = reshape(v,subn);
else
    v = reshape(v,model.n);
end
nlabs = parpool_size(); 

if exist('srcfreqmask','var')==0 || isempty(srcfreqmask)
    if isfield(params,'srcfreqmask')
        srcfreqmask = params.srcfreqmask;        
    else
        if ndims==2
            srcfreqmask = true(nsrc,nfreq);
        elseif nfreq==1
            srcfreqmask = true(nsrc,nfreq);
        else
            error('Need to specify srcfreqmask');
        end
    end
end
assert(size(srcfreqmask,1) == nsrc && size(srcfreqmask,2) == nfreq );
assert(islogical(srcfreqmask));    


nout = nargout;

% Serial version not recommended, you should do things in parallel for FWI
if nlabs == 0
    freqsxsy = srcfreqmask;
    if strcmp(src_est_mode,PDEopts.SRC_EST_NONE)==0 && strcmp(func,OBJ)
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

% Partition the requested (source, freq) pairs in a distributed manner
I = find(vec(srcfreqmask));
% Handle data distribution in 3D
if ndims==3    
    iF = find(any(srcfreqmask,1));
    assert(length(iF)==1,'For 3D, lower frequencies finish faster than higher frequencies, so loop your calls over frequency so your workers have balanced workloads');
    if ~isempty(Dobs) 
        assert(size(Dobs,2)==length(model.xsrc)*length(model.ysrc) || size(Dobs,2)==length(find(srcfreqmask(:,iF))));
        
        if size(Dobs,2)==length(model.xsrc)*length(model.ysrc) %Data is full freq slice, need to subsample + distribute, if necessary                                               
            if (~iscodistributed(Dobs)||~isdistributed(Dobs))                
                Dobs = Dobs(:,srcfreqmask(:,iF));
                Dobs = distributed(Dobs);                   
            else
                iS = distributed(srcfreqmask(:,iF));
                Dobs = Dobs(:,iS);
                spmd,
                    codist = getCodistributor(iS);
                    Dobs = redistribute(Dobs,codistributor1d(2,codist.Partition,size(Dobs)));
                end
            end            
        else %Data has already been subsampled, just need to distribute, possibly
            if (~iscodistributed(Dobs)||~isdistributed(Dobs))  
               Dobs = distributed(Dobs); 
            end
        end       
    end
end
if ~isempty(Dobs)
    spmd
        Dloc = getLocalPart(Dobs); 
    end
end

npdes = length(I);
spmd
    I = getLocalPart(codistributed(I));  
    freqsxsy = false(nsrc,nfreq);
    freqsxsy(I) = true;
    [~,iF] = ind2sub([nsrc,nfreq],I);
    ifreqs = min(iF):max(iF);
    if size(Q,3)>1 && (iscodistributed(Q)||isdistributed(Q))
        Q = getLocalPart(Q);
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
    warning('off','all');
    [~,host] = unix('hostname');
    switch func
      case OBJ       
        if strcmp(src_est_mode,PDEopts.SRC_EST_RECOMPUTE)
            input = {w(ifreqs)};
        elseif strcmp(src_est_mode,PDEopts.SRC_EST_NORECOMPUTE)
            input = {w(ifreqs),getLocalPart(U)};
        else        
            input = {};
        end                
        if nout == 1
            tic,f = PDEfunc(func,v,Q,input,Dloc,model,params,freqsxsy);T = toc;
            if disp_host                
                disp([host ' took ' num2str(T) 's']);
            end
        else
            
            tic,[f,g,h] = PDEfunc(func,v,Q,input,Dloc,model,params,freqsxsy);T = toc;
            if disp_host
            disp([host ' took ' num2str(T) 's']);
            end
            g = pSPOT.utils.global_sum(g,1);
            h = pSPOT.utils.global_sum(h,1);                 
        end            
        f = pSPOT.utils.global_sum(f,1); 
        dm = [];
      case JACOB_ADJ
        dm = PDEfunc(func,v,Q,getLocalPart(input),[],model,params,freqsxsy);
      case HESS
        dm = PDEfunc(func,v,Q,input,Dloc,model,params,freqsxsy);
      case HESS_GN
        dm = PDEfunc(func,v,Q,input,[],model,params,freqsxsy);        
      case {FORW_MODEL,JACOB_FORW}
        tic,D = PDEfunc(func, v, Q, input, [],model,params,freqsxsy);T = toc;
        if disp_host
            disp([host ' took ' num2str(T) 's']);
        end
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
                varargout = {f{1},g{1},h{1}};
            end
        end
    case {JACOB_ADJ,HESS,HESS_GN}
        varargout = {dm{1}};
        
    case {FORW_MODEL,JACOB_FORW}
        varargout = {D};
end
            

end

