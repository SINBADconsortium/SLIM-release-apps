function p = TraceNorm_primal_TJM(x, weights, params)

% Use:
%   TraceNorm_primal_TJM(x, weights, params)
%
% Input:   
%      x - unknow data 
%      weight - weight value       
%      params - parameter file

% Author: Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: January, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%----------------------------------------------------------------------------------------------------


p = 0.5 * norm(x.*weights)^2; 

end