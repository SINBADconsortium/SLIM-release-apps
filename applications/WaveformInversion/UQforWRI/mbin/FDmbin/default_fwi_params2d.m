function params = default_fwi_params2d(params)
% default_fwi_params2d - Appends default 2D FWI parameters to the
% input struct, if not specified. If your 2D FWI examples aren't
% working for some reason, you should think about specifying
% non-default options to your params struct.
% 
% Usage:
%   params = default_fwi_params2d(params);
%
% Input/output:
%   params - parameters struct passed on to
%   PDEfunc/PDEfunc_dist. The following defaults are specified if
%   they are not currently present.
% 
%   params.pdefunopts = PDEopts2D();
%   lsopts = LinSolveOpts();
%   lsopts.solver = LinSolveOpts.SOLVE_LU;
%   params.lsopts = lsopts;
%
% Curt Da Silva, 2016
    
    if isfield(params,'pdefunopts')==0
        params.pdefunopts = PDEopts2D();
    end
    if isfield(params,'lsopts')==0
        solve_opts = LinSolveOpts2D();
        params.lsopts = solve_opts;
    end
end