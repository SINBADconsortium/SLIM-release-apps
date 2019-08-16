function [D,J] = F3d(m,Q,model,params)
% 3D Frequency domain data modeling
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
    freqs = any(sf_mask,1);
    if nfreq == 1
        max_freq = 1;
    else                                  
        max_freq = find(freqs,1,'last');
        if isempty(max_freq), max_freq = length(freqs); end            
    end    
    min_freq = find(freqs,1,'first');
    if isfield(params,'subsample_model') && params.subsample_model
        params_loc = params;
        params_loc.pdefunopts = copy(params.pdefunopts);
        params_loc.lsopts = copy(params.lsopts);
        [dt,m] = subsample_model(m,model,params_loc.pdefunopts.helm_scheme,model.freq(max_freq));        
        params_loc.pdefunopts.helm_dt = dt;        
    else
        params_loc = params;
    end
    Ds = cell(length(freqs),1);
    ifreqs = find(freqs);
    for i=1:length(ifreqs)   
        sf_mask_f = false(nsrc,nfreq);
        sf_mask_f(:,ifreqs(i)) = sf_mask(:,ifreqs(i));
        Ds{i} = PDEfunc_dist(PDEopts.FORW_MODEL,m,Q,[],[],model,params_loc,sf_mask_f);
    end
    
    % Not the most efficient way to do this, but you only model once (yomo)
    if length(freqs)>1
        D = [];
        for i=1:length(freqs);
            D = [D,gather(Ds{i})];
        end
        D = distributed(D);
    else
        D = Ds{1};
    end
    if parpool_size()==0
        J = opDF3d(m,Q,model,params);
    else
        J = oppDF3d(m,Q,model,params);
    end
end