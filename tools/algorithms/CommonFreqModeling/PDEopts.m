classdef PDEopts < matlab.mixin.Copyable
%PDEOPTS - Options for PDEfunc
%   
%    Curt Da Silva, 2015
%
% Options
% 
% .numcompsrc           - number of sources to process at the same time (default: 4)
%
% .zero_boundary        - if true, zeros the gradient/hessian-vector product at the boundary nodes (default: false)
%
% .window_source_grad   - if true, zeros the gradient/hessian-vector product in a neighbourhood of the sources (default: false)
%
% .src_interp           - source interpolation method, one of 
%     PDEopts.SRC_INTERP_LIN  - linear interpolation
%     PDEopts.SRC_INTERP_SINC - sinc interpolation (default)
%
% .rec_interp           - receiver interpolation method, one of
%     PDEopts.REC_INTERP_LIN  - linear interpolation
%     PDEopts.REC_INTERP_SINC - sinc interpolation (default)
%
% .debug_mode           - if true, perform some basic debugging checks (default: false)  
%
% .helm_scheme          -  helmholtz discretization scheme, one of
%     PDEopts.HELM_OPERTO27   - 27 pt stencil based on Operto et al. (2007) (default 3D)
%     PDEopts.HELM_STD7       - 7 pt standard finite difference stencil
%     PDEopts.HELM2D_CHEN9P   - optimal 9pt stencil based on Chen et al. (2009) (default 2D)
%
% .helm_cut_pml         - if true (default), removes the contribution of the velocity model in the PML region, if false, uses the true adjoint of PML extension to return to the computational domain 
%
% .helm_free_surface    - if true, uses a free surface on the top of the computational domain, e.g., no pml (default: false)
%
% .helm_pml_max         - if specified, maximum number of pml points to add to the computational domain (default: inf)
%
% .helm_dt              - grid spacing for the computational domain (default: [], use model.d)
%
% .helm_pml             - either a scalar or a length 3 vector, specifying the number of pml points in the x-y-z coordinates (default: [], chosen by the modeling code)
%                          
% This class also defines the following variables for specifying the mode of PDEfunc / PDEfunc_dist
%   PDEopts.OBJ         - least squares objective/gradient
%          .FORW_MODEL  - forward modeling
%          .JACOB_FORW  - demigration-vector product
%          .JACOB_ADJ   - migration-vector product
%          .HESS_GN     - GN Hessian-vector product
%          .HESS        - full Hessian-vector product
    
    properties        
        numcompsrc, ...
        zero_boundary, ...
        window_source_grad, ...
        src_interp, ...
        rec_interp, ...
        debug_mode,...
        helm_scheme, ...
        helm_cut_pml,...
        helm_free_surface,...
        helm_pml_max,...
        helm_dt,...
        helm_pml, ...
        helm_mat_free,...
        helm_nthreads,...
        src_est_mode,...
        misfit_func,...
        offset_mask;        
    end
    
    properties (Constant)
        OBJ = 'obj';
        FIELD = 'field'
        JACOB_FORW = 'jacob_forw';
        JACOB_ADJ = 'jacob_adj';
        HESS = 'hess';
        HESS_GN = 'hess_gn';
        HESS_GN_DIAG = 'hess_gn_diag';
        HESS_DIAG_SHIN01 = 'hess_diag_shin01';
        HESS_DIAG_ENCODE = 'hess_diag_encode';
        HESS_NONE = 'hess_none';
        OP1 = 'op1';
        FORW_MODEL = 'forw_model';
        MODE_UNKNOWN = 'mode_unknown';
        SRC_INTERP_LIN = 'cubic';
        SRC_INTERP_SINC = 'sinc'; 
        REC_INTERP_LIN = 'cubic';
        REC_INTERP_SINC = 'sinc';
        HELM_OPERTO27 = 'operto27';
        HELM_STD7 = 'std7pts';
        HELM3D_OPERTO27 = 'operto27';
        HELM3D_STD7 = 'std7pts';
        HELM3D_CHEN27 = 'chen27';
        HELM3D_CHEN2012 = 'chen2012';
        
        HELM2D_JO9P = 'helm2d_9p';
        HELM2D_CHEN9P = 'helm2d_9popt';
        POISS2D_FV = 'poiss2d_finitevol';
        PML_EXT = 1;
        PML_NOEXT = 2;
        PML_NONE = 0;
        WAVELET_RICKER = 'wavelet_ricker';
        WAVELET_CUSTOM = 'wavelet_custom';
        SRC_EST_NONE = 'src_est_none';
        SRC_EST_RECOMPUTE = 'src_est_recompute';
        SRC_EST_NORECOMPUTE = 'src_est_norecompute';
        SRC_EST_LOC = 'src_est_loc';
    end
    
    methods
        function opts = PDEopts()            
            opts.numcompsrc         = 1;
            opts.zero_boundary      = false;
            opts.window_source_grad = 0;
            opts.src_interp         = PDEopts.SRC_INTERP_SINC;
            opts.rec_interp         = PDEopts.REC_INTERP_SINC;            
            opts.debug_mode         = false;
            opts.helm_scheme        = PDEopts.HELM3D_OPERTO27;
            opts.helm_cut_pml       = true;
            opts.helm_free_surface  = false;
            opts.helm_pml_max       = 10;
            opts.helm_dt            = [];
            opts.helm_pml           = 10;
            opts.helm_mat_free      = true;
            opts.helm_nthreads      = 5;
            opts.src_est_mode       = PDEopts.SRC_EST_NONE;
            opts.misfit_func        = @LS_misfit;
            opts.offset_mask        = [];
        end
    end
    
end

