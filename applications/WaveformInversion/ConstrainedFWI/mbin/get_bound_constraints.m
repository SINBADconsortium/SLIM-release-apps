function [LB,UB]=get_bound_constraints(m,comp_grid,constraint)
%sets up bound constraints, based on global min and max values, possibly specific bounds for the water layer and possible bound relative to a reference/start model
% Bas Peters

LB = ones(comp_grid.n)*constraint.v_min;
UB = ones(comp_grid.n)*constraint.v_max;

if isfield(constraint,'water_depth')==1
    water_bottom_index = floor(constraint.water_depth/comp_grid.d(1));
%     if water_bottom_index==0;
%         water_bottom_index=1;
%     end
    LB_water = ones(comp_grid.n).*constraint.water_min;
    UB_water = ones(comp_grid.n).*constraint.water_max;
    LB_water(water_bottom_index+1:end,:)=LB(water_bottom_index+1:end,:);
    UB_water(water_bottom_index+1:end,:)=UB(water_bottom_index+1:end,:);
    
    %now find tightest combination
    LB=max(LB(:),LB_water(:));
    UB=min(UB(:),UB_water(:));
    clear UB_water LB_water
end
if isfield(constraint,'relative_2_ini_plus')==1
    %now find tightest combination
    LB=max(LB(:),constraint.relative_2_ini_minus(:));
    UB=min(UB(:),constraint.relative_2_ini_plus(:));
end

LB = LB(:);
UB = UB(:);
end