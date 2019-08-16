function [D,J] = F3d(m,Q,model,params)
% 3D Frequency domain data modeling
% 
% Curt Da Silva, 2015
%
% Usage:
%   [D,J] = F3d(m,Q,model,params);
%
% Input:
%   m      - 3D model vector
%   Q      - source weight matrix
%   model  - model parameters struct
%   params - performance parameters struct
%
% Output:
%   D      - forward modelled data
%   J      - pSPOT demigration/migration operator
    nsrc = size(Q,2); nfreq = length(model.freq);    
    if ~isfield(params,'srcfreqmask'), params.srcfreqmask = true(nsrc,nfreq); end
    sf_mask = params.srcfreqmask;

    if isfield(params,'subsample_model') && params.subsample_model
        freqs = any(sf_mask,1);
        if nfreq == 1
            max_freq = 1;
        else                                  
            max_freq = find(freqs,1,'last');
            if isempty(max_freq), max_freq = length(freqs); end            
        end
        params_loc = params;
        params_loc.pdefunopts = copy(params.pdefunopts);
        params_loc.lsopts = copy(params.lsopts);
        [dt,m] = subsample_model(m,model,params_loc.pdefunopts.helm_scheme,model.freq(max_freq));        
        params_loc.pdefunopts.helm_dt = dt;        
    else
        params_loc = params;
    end
    D = PDEfunc_dist(PDEopts.FORW_MODEL,m,Q,[],[],model,params_loc);
    if parpool_size()==0
        J = opDF3d(m,Q,model,params);
    else
        J = oppDF3d(m,Q,model,params);
    end
end