function model_bg = smooth_ND_model(model_true,k,smooth_slowness_flag)
% Syntax
% vel_bg = smooth_ND_model(vel_true,k,smooth_slowness_flag)
%
% Description
% smooth the true model to obtain a backgroud smooth model.
%
% Input list:
% model_true: true model
% k: the number of times you apply the moving average filter [0.25
% 5 0.25]
% smooth_slowness_flag:
% 0, smooth velocity
% 1, smooth slowness
% 2, smooth slowness square
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

if nargin < 3
    % by default smooth slowness square
    smooth_slowness_flag = 2;
end

size_model = size(model_true);
N = length(size_model);
op_smooth_all = cell(N,1);

for count = 1:N
    op_smooth_all{N+1-count} = opSmooth(size_model(count),k);
end

op_smooth_all = opKron(op_smooth_all{:});
switch smooth_slowness_flag
    case 0
        model_bg = op_smooth_all*model_true(:);
    case 1
        model_bg = 1./(op_smooth_all*vec(1./model_true));
    case 2
        model_bg = 1./sqrt(op_smooth_all*vec(1./(model_true.^2)));
end
model_bg = reshape(model_bg,size_model);