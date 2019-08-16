function d = TraceNorm_dual(x, weights, params)

%-----------------------------------------------------------------------
% Dual of trace norm is operator norm, i.e., the maximum singular value
%
% Use:
%   TraceNorm_dual_macq(x, weights, params)

% Author: Rajiv Kumar and Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%-----------------------------------------------------------------------

fps1 = reshape(x(1:params.mhnumr*params.mhnumc),params.mhnumr,params.mhnumc);
fps2 = reshape(x(params.mhnumr*params.mhnumc+1:end),params.mhnumr,params.mhnumc);
E = [fps1;fps2];
d = svds(gather(E),1);    

end % function end

