%-------------------------------------------------------------------------------
% Discard the receivers which are closer than hmin to the source.
% This returns a binary diagonal sparse matrix. Performing a matrix-vector
% product between this matrix and the data vector removes the said receivers.
% 
% Setting hmin to -1 (or any negative number) will return the identity.
% 
% Use: 
%   Wh = offset_mask(acq,hmin,l)
% Input:
%  acq        - Contains information about the acquisition and its geometry. 
%               Read further for more details on its content.
%  hmin       - Minimum offset to be considered
%  l          - The output sparse diagonal matrix Wh contains the proper
%               weighting for the l-th source.
%
% Output:
%   Wh        - Sparse diagonal matrix containing the proper weighting.
% 
% Structure of "acq" (input):
% ---------------------------------------
% *Must* contain:
% 
%  acq.{zsrc,xsrc,ysrc} - vectors describing source array
%  acq.{zrec,xrec,yrec} - vectors describing receiver array.
%
% Author: Tristan van Leween
%         Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: September, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%-------------------------------------------------------------------------------
function Wh = offset_mask(acq,hmin,l)
   
   % Retrieve the 3d index of the l-th source:
   [ix,iy,iz]=ind2sub([length(acq.xsrc) length(acq.ysrc) length(acq.zsrc)],l);
   
   % Obtain distance between the said source and the receivers
   [hx,hy,hz]=ndgrid(acq.xrec-acq.xsrc(ix),acq.yrec-acq.ysrc(iy),acq.zrec-acq.zsrc(iz));

   % Check if it is within the oftset range.
   Wh = sqrt((length(acq.yrec)>1)*hy.^2 ... 
           + (length(acq.xrec)>1)*hx.^2 ...
           + (length(acq.zrec)>1)*hz.^2) > hmin;

   Wh=spdiags(Wh(:), length(Wh(:)),length(Wh(:)));
end 
