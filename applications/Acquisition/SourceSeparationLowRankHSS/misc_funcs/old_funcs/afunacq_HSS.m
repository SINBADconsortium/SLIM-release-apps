function [output] = afunacq_HSS(x, MH, params, mode)

%------------------------------------------------------------------
% Use:
%   afunacq_HSS(x, MH, params, mode)
%
% Input:   
%          x - input data vector
%         MH - midpoint-offset operator
%     params - parameters set in the SourceSep_params.m and SourceSep_HSS.m scripts
%       mode - [1 : forward], [-1 : adjoint]

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------

if mode == 1 
   output = opBlockDiag(MH',MH')*vec(x);  % from MH to SR
   output = shiftzs_HSS(output,params,1);
else
   output = shiftzs_HSS(x,params,-1);
   output = opBlockDiag(MH,MH)*vec(output);  % from SR to MH with shift
end

end % function end

