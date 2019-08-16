function [AT, Dden] = selectCoeffSupp_denoise(C, xL1, perc, D, nmfactor)

%------------------------------------------------------------------------------
% selectCoeffSupp_denoise selects the support of the significant curvelet coefficients that represent the desired seismic signal and denoises
%
% Use:
%   [AT, Dden] = selectCoeffSupp_denoise(C, xL1, perc, D, nmfactor)
%
% Input: 
%            C - 2D curvelet operator
%          xL1 - synthesis curvelet coefficients (output of getCurveletCoeff.m)
%         perc - percentage of energy captured by the significant curvelet coefficients
%            D - input (noisy) data
%     nmfactor - data normalization factor  
%
% Output:
%        AT - restricted curvelet operator (restricted set of curvelets that we identified to significantly represent our desired seismic signal)
%      Dden - denoised data
      
% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: August, 2012

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------------------

A = C';
[Cx idx] = sort(abs(xL1), 'descend');
ratCx = sqrt(cumsum(Cx.^2))/norm(Cx);
k = find(ratCx >= perc);
AT = A(:, idx(1:k(1)));
xden = pinv(AT)*D(:);
dim = size(D);
Dden = reshape(AT*xden, dim); 
Dden = Dden*nmfactor;

end  % function

