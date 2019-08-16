function [obj,vout,comp_grid] = misfit_setup(v0,Q,D,model,params)
% MISFIT_SETUP - Generates a least-squares FWI or WRI objective for use in an optimization algorithm.
% Handles data parallelization along frequencies and sources automatically.  
%
% Usage:
%   func = misfit_setup(v0,Q,D,model,params);
%
% Input:
%   Q     - source matrix of size nsrc x number_of_sources_used
%   D     - vectorized, distributed data of size nrec x nsrc x nfreq
%   model - struct with model parameters
%   params 
%   
%       .wri             - if true (default), computes WRI objective,
%                          if false, computes standard FWI objective
%       .lambda          - tradeoff parameter (scalar) between PDE and data
%                          misfit (for WRI only)
%       .dist_mode       - 'freq'    : distributes over frequencies
%                          'srcfreq' : distributes over sources and frequencies
%                            
%       .srcfreqmask     - nsrc x nfreq binary mask
%                          1 at (i,j)   - compute with source i, freq j
%                          otherwise    - ignore source i, freq j
%
%       .hessian         - 'sparse' : WRI only, sparse Hessian (default for WRI)
%                          'gn'     : Gauss-Newton Hessian (default for FWI)
%                          'full'   : Full Hessian 
%
%       .subsample_model - if true, coarsen the model to the number of points per wavelength specified by the stencil (default: false)
%       
%       .batch_mode      - if true, function handle accepts source/freq indices to compute with at function call time, see below (default: false)
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
%   if params.batch_mode is true, evaluate the objective as
%     [fk,gk,hk] = func(m,I);
%   where I is a vector of integers with values anywhere from 1 to length(find(params.srcfreqmask)) (i.e., number of sources and frequencies used)
%
%  Author: Curt Da Silva, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
    mode_2d = length(model.n)==2 || model.n(3)==1;
    if mode_2d
        if exist('params','var')==0, params = struct; end
        params = default_fwi_params2d(params);
        model = fwi2d_model_compatibility(model);        
        if length(model.xrec) > 1
            nrec = length(model.xrec);
        else
            nrec = length(model.zrec);
        end
    else
        if exist('params','var')==0, error('Need params struct'); end
        params.wri = false;        
        nrec = length(model.xrec)*length(model.yrec);
    end    
    
    nsrc = size(Q,2);  
    nfreq = length(model.freq);

    if exist('params','var')==0, error('Need to specify params'); end
    if ~isfield(params,'hessian'), params.hessian = PDEopts.HESS_GN; end
    if ~isfield(params,'srcfreqmask'), params.srcfreqmask = ones(nsrc,nfreq); 
    else assert(norm(size(params.srcfreqmask)-[nsrc,nfreq])==0, 'src freq mask must have size nsrc x nfreq'); end
    if isfield(params,'subsample_model'),sub_model = params.subsample_model; else sub_model = false; end
    if isfield(params,'batch_mode'),batch_mode = params.batch_mode; else batch_mode = false; end    
            
    sf_mask = params.srcfreqmask;
    
    freqs = any(sf_mask,1); 
    nfreq = length(find(freqs));
    
    if nfreq == 1
        min_freq = 1; max_freq = 1;
    else
        % All of the frequencies before this index are 0
        min_freq = find(freqs,1,'first');
        if isempty(min_freq), min_freq = 1; end
        
        % All of the frequencies after this index are 0
        max_freq = find(freqs,1,'last');
        if isempty(max_freq), max_freq = length(freqs); end
    end
    params_loc = params;
    params_loc.pdefunopts = copy(params.pdefunopts);
    params_loc.lsopts = copy(params.lsopts);
    
    comp_grid = struct;
    
    % Subsample model, if requested
    if sub_model           
        [dt,vout,subn,to_coarse,to_fine] = subsample_model(v0,model,params_loc.pdefunopts.helm_scheme,model.freq(max_freq));
        params_loc.pdefunopts.helm_dt = dt;
    else
        to_fine = opDirac(numel(v0));
        to_coarse = opDirac(numel(v0));
        dt = model.d; subn = model.n; 
        vout = v0;
    end
    
    comp_grid.dt = dt; comp_grid.nt = subn; comp_grid.ot = model.o;
    comp_grid.to_coarse = to_coarse;
    comp_grid.to_fine = to_fine;
    
    % Parse data, either we get a function handle pointing to the paths of the data or the full data itself
    if isa(D,'function_handle')
        fslice_name_func = D;
        I = distributed(find(sf_mask));
        num_pts = length(I);
        spmd
            codist = getCodistributor(I);
            I = getLocalPart(I);
            [idx_src,idx_freq] = ind2sub([nsrc,nfreq],I);
            freqs_loc = unique(idx_freq);
            D1 = zeros(nrec,length(I));
            ns = 1;
            for i=1:length(freqs_loc)
                fi = freqs_loc(i);
                filepath = fslice_name_func(fi);
                Dloc = load(filepath);
                fields = fieldnames(Dloc);
                if length(fields) > 1, error('Too many variables in file'); end
                Dloc = Dloc.(fields{1});
                assert(size(Dloc,1)==nrec && size(Dloc,2)==nsrc,['Unexpected freq slice size, expected ' num2str(nrec) ' x ' num2str(nsrc) ' but got ' num2str(size(Dloc,1)) ' x ' num2str(size(Dloc,2))]);
                iS_f = idx_src(idx_freq==fi);
                D1(:,ns:ns+length(iS_f)-1) = Dloc(:,iS_f);
                ns = ns+length(iS_f);
            end
            codistr = codistributor1d(2,codist.Partition,[nrec,num_pts]);
            Dobs = codistributed.build(D1,codistr,'noCommunication');
        end    
    elseif isdistributed(D) || iscodistributed(D)
        if numel(D)==nrec*nsrc*length(model.freq)
            nf = length(model.freq);
        elseif numel(D)==nrec*nsrc*nfreq
            nf = nfreq;
        else
            error('Unrecognized size for D');
        end
        
        Dobs = reshape(D,[nrec,nsrc*nf]);
        % If we do any subsampling, subsample data according to specified src-freq mask + redistribute it 
        if sum(vec(sf_mask)) ~= numel(sf_mask)
            %If it's a single frequency slice requested + passed in, don't do anything
            if min_freq~=max_freq        
                Dobs = Dobs(:,vec(sf_mask));
                spmd,
                    Dobs = redistribute(Dobs,codistributor1d(2,codistributor1d.unsetPartition,size(Dobs)));
                end 
            end
        end
    else
        error('Unrecognized data type D -- must be a distributed array ' );
    end    
    
    % Output objective   
    if batch_mode
        obj = @(m,I) misfit_func(reshape(m,comp_grid.nt),Q,Dobs,model,params_loc,I);
    else
        obj = @(m) misfit_func(reshape(m,comp_grid.nt),Q,Dobs,model,params_loc);
    end    
end

function [f,g,h] = misfit_func(m,Q,D,model,params,I)
    if exist('I','var')==0
        I = 1:length(find(vec(params.srcfreqmask)));
    end
    % Source/freq subsampling can only be done within the current
    % (fixed) set of source/freq indices
    Ifixed = find(vec(params.srcfreqmask));
    assert(max(I) <= length(Ifixed),'Requested index set out of range of current sources/frequencies');
    sfmask = false(size(params.srcfreqmask));
    sfmask(Ifixed(I)) = true; 
    params.srcfreqmask = sfmask;
    
    if length(I) < size(D,2)
        If = distributed(I);
        D1 = D(:,If);
        spmd,
            codist = getCodistributor(If);
            D1 = redistribute(D1,codistributor1d(2,codist.Partition,size(D1)));
        end
    else
        D1 = D;
    end    
        
    c = 1/length(I);
    if params.wri       
        [f,g,h] = misfit_pen(m,Q,D1,model,params);
        f = c*f; g = c*g; h = c*h;
        switch params.hessian
          case 'gn'
            h = c*oppHGNpen(m,Q,D1,model,params);
          case 'full'
            h = c*oppHpen(m,Q,D1,model,params);
        end
    else    
        if nargout== 1
            f = PDEfunc_dist(PDEopts.OBJ,m,Q,[],D1,model,params);
        else
            [f,g] = PDEfunc_dist(PDEopts.OBJ,m,Q,[],D1,model,params);
            g = c*g;        
        end
        f = c*f;
        if nargout >= 3
            switch params.hessian
              case PDEopts.HESS_GN
                h = c*oppHGN3d(m,Q,model,params);
              case PDEopts.HESS
                h = c*oppH3d(m,Q,D1,model,params);
            end        
        end
    end
    
end