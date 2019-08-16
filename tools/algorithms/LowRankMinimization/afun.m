function [y] = afun(x,Ind)
% Insert zeros back in the interpolated matrix to compute the residual.
% 
% use:
%   [y] = afun(x,Ind)
%
% input:
% x         - input data
% Ind       - Locations of zeros corresponds to the missing data in observed data.
% output:
%  y        - output with zeros at pre-defined locations based on the
%             observed data.
% Author: Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
 x = vec(x);
 x(Ind)=0;
 y = x;
 end