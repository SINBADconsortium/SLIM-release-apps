function [TD_op]=get_TD_operator(comp_grid,type)

%if isfield(constraint,'card_operator'); type=constraint.card_operator; end
%if isfield(constraint,'TD_operator');   type=constraint.TD_operator;   end

h1 = comp_grid.d(1);
h2 = comp_grid.d(2);
n1 = comp_grid.n(1);
n2 = comp_grid.n(2);

if strcmp(type,'TV') %% TV operator
    [TD_op,~] = getDiscreteGrad(n1,n2,h1,h2,ones(n1,1));
elseif strcmp(type,'TV2')%% vertical derivative operator
    [TD_op,~] = getDiscreteGrad(n1,n2,h1,h2,ones(n1,1));
    TD_op=TD_op'*TD_op;
elseif strcmp(type,'D_vert')%% vertical derivative operator
    [S,~] = getDiscreteGrad(n1,n2,h1,h2,ones(n1,1));
    s = (n1-1)*(n2);
    TD_op=S(1:s,:);
    
elseif strcmp(type,'D_lat')%% lateral derivative operator
    [S,~] = getDiscreteGrad(n1,n2,h1,h2,ones(n1,1));
    s = (n1-1)*(n2);
    TD_op=S(s+1:end,:);
    
elseif strcmp(type,'curvelet')%% curvelet
    TD_op=opCurvelet(n1,n2);
elseif strcmp(type,'DFT2D')%% DFT 2D
    TD_op=opDFT2(n1,n2);
elseif strcmp(type,'DFT2D_plus_curvelet')%% DFT 2D
    TD_op=[opDFT2(n1,n2);opCurvelet(n1,n2)];
else
    error('provided an unknown transform domain operator. check function [TD_op]=get_TD_operator(comp_grid,type) for options')
end
end