function [TD_op]=get_TD_operator(m,model_PDE,constraint)

if isfield(constraint,'card_operator'); type=constraint.card_operator; end
if isfield(constraint,'TD_operator');   type=constraint.TD_operator;   end
h1=model_PDE.d(1);
h2=model_PDE.d(2);

if strcmp(type,'TV') %% TV operator
    [TD_op,~] = getDiscreteGrad(model_PDE.n(1),model_PDE.n(2),h1,h2,ones(size(m,1),1));
    
elseif strcmp(type,'D_vert')%% vertical derivative operator
    [S,~] = getDiscreteGrad(model_PDE.n(1),model_PDE.n(2),h1,h2,ones(size(m,1),1));
    s = (model_PDE.n(1)-1)*(model_PDE.n(2));
    TD_op=S(1:s,:);
    
elseif strcmp(type,'D_lat')%% lateral derivative operator
    [S,~] = getDiscreteGrad(model_PDE.n(1),model_PDE.n(2),h1,h2,ones(size(m,1),1));
    s = (model_PDE.n(1)-1)*(model_PDE.n(2));
    TD_op=S(s+1:end,:);
    
elseif strcmp(type,'curvelet')%% curvelet
    TD_op=opCurvelet(model_PDE.n(1),model_PDE.n(2));
end
end