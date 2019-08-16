function update = Smoothedge(up,depth,mode)
% Usually in the model update, there are a lot of informations in the water layer.
% which should not be there. This function gits rid of the information in the water
% layer.
% 
% use: update = Smoothedge(up,depth,mode)
% 
% Input:
% 	up: image
% 	depth: approximated bottom of the water layer
% 	mode: if 1, smooth the edge of bottom; if 2, hard cut. 
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

mask = ones(size(up,1),1);
mask(1:depth) = 0;
mask1 = smooth(smooth(mask));

if mode == 1
	% smooth mode
	update = repmat(mask1,1,size(up,2)).*up;
elseif mode == 2
	% hard mode
	update = repmat(mask,1,size(up,2)).*up;
end