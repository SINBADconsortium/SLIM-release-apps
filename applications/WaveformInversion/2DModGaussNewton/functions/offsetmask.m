function [index mask] = offsetmask(nx,Grid_length,sx,minoffst,maxoffst);
%  This is a function which can create a mask for the wavefield to
%  get rid of near offset and far offset data.
%
% Input:
% 	nx: number of lateral grid points of the model
% 	Grid_length: lateral mesh distance
% 	sx: source point position
% 	minoffst: minimul offset
% 	maxoffst: maximul offset
%
% 
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


minp = round(minoffst./Grid_length) - 1;
maxp = round(maxoffst./Grid_length) - 1;
ns   = length(sx);
mask = zeros(nx,ns);

for m = 1:ns
	index =  max([sx(m)-maxp,1]) : min([sx(m) + maxp,nx]);
	mask(index,m) = 1;
	index =  max([sx(m)-minp,1]) : min([sx(m) + minp,nx]);
	mask(index,m) = 0;
end


index = find(mask~=0);



