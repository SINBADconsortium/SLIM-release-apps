function [dt,vc,subn,to_coarse,to_fine] = subsample_model(v,model,helm_scheme,max_freq)
% subsample_model - Subsample a velocity model for a helmholtz equation according to the maximum allowed spacing of the helmholtz scheme
%
% Curt Da Silva, 2016
%
% Usage:
%   [dt,vc,subn,to_coarse,to_fine] = subsample_model(v,unit,helm_scheme,max_freq);
%
% Input:
%   v           - velocity model
%   unit        - units of v ('m/s','s2/m2', or 's2/km2')
%   helm_scheme - helmholtz scheme
%   max_freq    - maximum frequency this velocity model will be computed on
% 
% Output:
%   dt          - coarsened model spacing
%   vc          - coarsened model
%   subn        - coarsened model dimensions
%   to_coarse   - fine grid to coarse grid operator
%   to_fine     - coarse grid to fine grid operator
    
    nmin = 50;
    dt = helm_stable_d(v,model.unit,helm_scheme,max_freq);
    if max(dt./model.d) < 1, error('Current model spacing potentially unstable for wave simulation, aborting'); end
    
    subn = ceil( model.n .* model.d ./ dt);
    J = find(subn <= nmin);
    if ~isempty(J)
        dt(J) = dt(J)./(nmin./subn(J));
    end
    [to_coarse,to_fine,subn] = fine2coarse(model.n,model.d,dt);        
    
    vc = to_coarse*vec(v);
end