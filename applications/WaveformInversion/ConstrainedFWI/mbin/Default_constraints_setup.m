function [ constraint,options ] = Default_constraints_setup()
%Dykstra options
constraint.options_dyk.maxIt=20;
constraint.options_dyk.minIt=3;
constraint.options_dyk.evol_rel_tol=1e-5;
constraint.options_dyk.tol=1e-2;
constraint.options_dyk.feas_tol=1e-1;
constraint.options_dyk.log_feas_error=0;
constraint.options_dyk.log_vec=0;

constraint.ADMM_options.maxit=800;
constraint.ADMM_options.evol_rel_tol=1e-6;
constraint.ADMM_options.rho=1e2;
constraint.ADMM_options.adjust_rho=1;
constraint.ADMM_options.test_feasibility=0;
constraint.ADMM_options.feas_tol=1e-2;

%% optimization options

options.maxIter = 15;
options.memory  = 5;
options.testOpt = 0;
options.verbose = 0;
options.interp=0;

%Total Variation
constraint.TV_initial    = 4;
constraint.TV_increment  = 2;
constraint.TV_proj       = 'admm';
constraint.TV_units      = 's^2/km^2';%'s^2/m^2';
% Box constraint
constraint.v_max   = 6; %upper bound
constraint.v_min   = 1.5; %lower bound
constraint.v_max_plus   = 1.5; %upper bound is current highest velocity + params.v_max_plus
constraint.v_min_minus  = .05; %lower bound is current lowest velocity -  params.v_min_minus
% minimum smoothness    full none aspect
constraint.smoothpars    = [4, 5, 3];%(0.8)

%% Se unused constraints to 0 
constraint.use_rank=0;
constraint.use_nuclear=0;


end

