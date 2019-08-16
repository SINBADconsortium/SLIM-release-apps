function [x] = TraceNorm_project_TJM(x,weights, B,params)
	
% Use:
%   TraceNorm_project_TJM(x,weights, B,params)
%
% Input:   
%      x - unknow data 
%      weight - weight value       
%      params - parameter file
%      B - tau value
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

c=sqrt(B/(0.5*norm(x)^2));
x = min(1,c)*x;

end