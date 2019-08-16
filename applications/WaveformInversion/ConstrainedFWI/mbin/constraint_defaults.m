function [constraint]=constraint_defaults(constraint)
%this functions loads all required fields into the 'constraint' structure
%if they are not provided

%select constraints defaults
    constraint.use_bounds       = 0;
    constraint.use_min_smooth   = 0;
    constraint.use_nuclear      = 0;
    constraint.use_rank         = 0;
    constraint.use_TV           = 0;
    constraint.use_cardinality  = 0;
    constraint.use_cardinality_2  = 0;
    constraint.use_TD_bounds    = 0;
    constraint.use_cvxbody      = 0;
    constraint.use_subspace     = 0;
    
end