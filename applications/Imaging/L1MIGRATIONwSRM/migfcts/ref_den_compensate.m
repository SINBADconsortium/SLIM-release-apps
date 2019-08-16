function model = ref_den_compensate(model, window_start, window_end, compensate_factor, taper_ratio)
% Syntax
% model_pert = ref_den_compensate(model_pert, window_start, window_end, compensate_factor)
%
% Description
% Increase ocean bottom reflectivity for multiple generation.
% If you use a modeling package assuming constant density, acoustic
% impedance introduced by density variation will not be taken into account.
% This code is to compensate the reflectivity caused by density variation
% by amplify the reflectivity by 'compensate factor'.
%
% Input list:
% model: model perturbation
% window_start: start position of ocean bottom reflector
% window_end: end of ocean bottom reflector
% compensate_factor: how much you want to amplity the reflector
% taper_ratio: a taper is applied to the ocean-bottom selection
% window, use taper_ratio to define how much taper you want to
% apply to the window.
%
% Output list:
% model_bg: the smooth model
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012

if nargin < 5
    % can be treated as alpha in a Tukey window function
    taper_ratio = 0.2;
end
size_model = size(model);

% decide window
model_ob_part = model(window_start:window_end,:);
length_ob = window_end - window_start +1;
Taper_ob = opKron(opDirac(size_model(2)),opWindow(length_ob,'Tukey',taper_ratio));

% compensate
model_ob_part = (compensate_factor-1)*Taper_ob*model_ob_part(:);
model_ob_part = reshape(model_ob_part,[length_ob,size_model(2)]);
model(window_start:window_end,:) = model(window_start:window_end,:) + model_ob_part;