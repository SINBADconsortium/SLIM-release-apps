function disp_verb(content, flag)
% Syntax:
% disp_verb(content, flag)
%
% Description:
% display with verbose flag
%
% Input list:
% content: the message to display
% flag: verbose flag
%
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: May/07/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if not(exist('flag','var'))
    flag = 1;
end
if flag
    disp(content);
end