function [constraint]=constraint_defaults(constraint)
%this functions loads all required fields into the 'constraint' structure
%if they are not provided
% Bas Peters, 2016
projection_accuracy = 'medium';
if exist('constraint') && isfield(constraint,'projection_accuracy'); projection_accuracy = constraint.projection_accuracy; end;


%select constraints defaults
constraint.use_bounds       = 0;
constraint.Minkowski1.use_bounds = 0;
constraint.Minkowski2.use_bounds = 0;
constraint.use_min_smooth   = 0;
constraint.use_nuclear      = 0;
constraint.use_rank         = 0;
constraint.use_TV           = 0;
constraint.use_cardinality  = 0;
constraint.use_cardinality_2= 0;
constraint.use_TD_bounds    = 0;
constraint.use_TD_l1        = 0;
constraint.use_TD_l2        = 0;
constraint.use_TD_l2_2      = 0;
constraint.use_cvxbody      = 0;
constraint.use_subspace     = 0;
constraint.use_quad_ineq    = 0;
constraint.use_poscard      = 0;
constraint.use_col_card     = 0;
constraint.use_hinge        = 0;

switch projection_accuracy
    case 'medium'
        %Dykstra options
        constraint.options_dyk.maxIt=8;
        constraint.options_dyk.minIt=1;
        constraint.options_dyk.evol_rel_tol=1e-4;
        constraint.options_dyk.tol=1e-2;
        constraint.options_dyk.feas_tol=1e-1;
        constraint.options_dyk.log_feas_error=0;
        constraint.options_dyk.log_vec=0;
        
        %ADMM options
        constraint.ADMM_options.maxit=500;
        constraint.ADMM_options.evol_rel_tol=2e-5;
        constraint.ADMM_options.rho=1e2;
        constraint.ADMM_options.adjust_rho=1;
        constraint.ADMM_options.test_feasibility=0;
        constraint.ADMM_options.feas_tol=1e-2;
        
    case 'budget'
        %Dykstra options
        constraint.options_dyk.maxIt=5;
        constraint.options_dyk.minIt=1;
        constraint.options_dyk.evol_rel_tol=1e-3;
        constraint.options_dyk.tol=2e-2;
        constraint.options_dyk.feas_tol=1e-1;
        constraint.options_dyk.log_feas_error=0;
        constraint.options_dyk.log_vec=0;
        
        %ADMM options
        constraint.ADMM_options.maxit=400;
        constraint.ADMM_options.evol_rel_tol=2e-4;
        constraint.ADMM_options.rho=1e2;
        constraint.ADMM_options.feas_tol=1e-2;
    case 'accurate'
        %Dykstra options
        constraint.options_dyk.maxIt=15;
        constraint.options_dyk.minIt=1;
        constraint.options_dyk.evol_rel_tol=2e-5;
        constraint.options_dyk.tol=1e-3;
        constraint.options_dyk.feas_tol=5e-2;
        constraint.options_dyk.log_feas_error=0;
        constraint.options_dyk.log_vec=0;
        
        %ADMM options
        constraint.ADMM_options.maxit=1500;%750;
        constraint.ADMM_options.evol_rel_tol=1e-6;%5e-6;
        constraint.ADMM_options.rho=1e1;
        constraint.ADMM_options.adjust_rho=1;
        constraint.ADMM_options.test_feasibility=0;
        constraint.ADMM_options.feas_tol=1e-3;
end

end