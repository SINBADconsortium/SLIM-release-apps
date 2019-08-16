function op = opRealRestriction(n)
% Syntax:
% op = opRealRestriction(n)
%
% Description:
% Remove the imaginary part of a vector
%
% Input list:
% n: length of vector
%
% Output list:
% op: the operator
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

funhandle = @(x,mode) opReal_intrnl(n,x,mode);
op=opFunction(n,n,funhandle);


function y = opReal_intrnl(n,x,mode)

if (mode == 0)
   y = {n,n,[0,0,0,0],{'RealRestriction'}};
else
   y = real(x);
end
