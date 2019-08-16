classdef PDEopts2D < PDEopts
%PDEOPTS - Default options for 2D FWI
%   
%    Curt Da Silva, 2016
%
% Options
% 
%
%
% .src_interp           - source interpolation method, one of 
%     PDEopts.SRC_INTERP_SINC - sinc interpolation (default)
%
% .rec_interp           - receiver interpolation method, one of
%     PDEopts.REC_INTERP_SINC - sinc interpolation (default)
%
% .debug_mode           - if true, perform some basic debugging checks (false)  
%
% .helm_scheme          -  helmholtz discretization scheme, one of
%     PDEopts.HELM2D_CHEN9P    - optimal 9pt stencil based on Chen et al. (2009)
%
% .helm_cut_pml         - if true (default), removes the contribution of the velocity model in the PML region (true)
%
% .helm_free_surface    - if true, uses a free surface on the top of the computational domain, e.g., no pml (false)
%
% .helm_pml_max         - if specified, maximum number of pml points to add to the computational domain (default: inf)
%
% .helm_dt              - grid spacing for the computational domain (default: [], use model.d)
%
% .helm_pml             - either a scalar specifying the number of pml points in the x-y-z coordinates (default: [], chosen by the modeling code)
%                          
    
    properties (Constant)
        default_numcompsrc = inf;
        default_zero_boundary = false;
        default_window_source_grad = false;
        default_source_est = PDEopts.SRC_EST_NONE;        
        default_src_interp = PDEopts.SRC_INTERP_SINC;
        default_rec_interp = PDEopts.REC_INTERP_SINC;
        default_debug_mode = false;
        default_helm_scheme = PDEopts.HELM2D_CHEN9P;
        default_helm_cut_pml = true;
        default_helm_free_surface = false;
        default_helm_pml_max = inf;
        default_helm_pml = 20;
        default_helm_mat_free = false;
    end
    
    methods
        function opts = PDEopts2D()            
            opts.numcompsrc         = PDEopts2D.default_numcompsrc;
            opts.zero_boundary      = PDEopts2D.default_zero_boundary;
            opts.src_est_mode       = PDEopts2D.default_source_est;
            opts.window_source_grad = PDEopts2D.default_window_source_grad;
            opts.src_interp         = PDEopts2D.default_src_interp;
            opts.rec_interp         = PDEopts2D.default_rec_interp;
            opts.debug_mode         = PDEopts2D.default_debug_mode;
            opts.helm_scheme        = PDEopts2D.default_helm_scheme;
            opts.helm_cut_pml       = PDEopts2D.default_helm_cut_pml;
            opts.helm_free_surface  = PDEopts2D.default_helm_free_surface;
            opts.helm_pml_max       = PDEopts2D.default_helm_pml_max;
            opts.helm_dt            = [];
            opts.helm_pml           = PDEopts2D.default_helm_pml;
            opts.helm_mat_free      = PDEopts2D.default_helm_mat_free;
        end
    end
    
end

