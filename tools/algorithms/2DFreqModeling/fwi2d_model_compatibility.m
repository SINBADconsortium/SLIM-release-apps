function model = fwi2d_model_compatibility(model)
% Conform older 2D FWI model structs to the latest requirements
%
% Curt Da Silva, 2016
%
    if isfield(model,'unit')==0, model.unit = 's2/km2'; end
    if length(model.o)==2, model.o = [model.o(:)',0]; end
    if length(model.d)==2, model.d = [model.d(:)',1]; end
    if length(model.n)==2, model.n = [model.n(:)',1]; end
end