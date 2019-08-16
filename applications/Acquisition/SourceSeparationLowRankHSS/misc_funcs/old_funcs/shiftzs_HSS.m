function [output] = shiftzs_HSS(x, params, mode)

%------------------------------------------------------------------
% shiftzs_HSS applies the time delay(s) and vice-versa.
%
% Use:
%   shiftzs_HSS(x, params, mode)
%
% Input:   
%          x - input data vector 
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

time = 2*params.depth/params.watervel;
shiftforward = exp(-1i*params.omega*(time + params.tdelay_sub));
shiftadjoint = exp(1i*params.omega*(time + params.tdelay_sub));

if mode == 1
   x1 = x(1:params.srnumr*params.srnumc);
   x2 = x(params.srnumr*params.srnumc+1:end);
   S1 = x1;
   S2 = shiftforward(:).*x2;
   output = S1 + S2; 
else
   S1 = x;
   S2 = shiftadjoint(:).*x;      
   output = [S1;S2];
end

end % function end

