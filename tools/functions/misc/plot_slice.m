% Plots a slice of a 3D cube.
% 
% Use:
%   plot_slice(V,s,d,type)
%
% Input:
%  V     - the thing you're trying to plot
%  s     - the index of the slice you're trying to plot, 
%  d     - the direction you're trying to slice
%  type  - its either "img" for image or "surf" for surface. If not given,
%          it will plot an image.
%          
% 
% Author: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: May, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%-----------------------------------------------------------------------------
function plot_slice(V,s,d,type)
if ~exist('type','var')
   type='img';
end

[nx, ny, nz] = size(V);

if strcmp(type,'img')
   if strcmp(d,'z')
      imagesc(reshape(V(:,:,s),nx,ny)');
   elseif strcmp(d,'y')
      imagesc(reshape(V(:,s,:),nx,nz)');
   elseif strcmp(d,'x')
      imagesc(reshape(V(s,:,:),ny,nz)');
   end
elseif strcmp(type,'surf')
   if strcmp(d,'z')
      surf(reshape(V(:,:,s),nx,ny)');
   elseif strcmp(d,'y')
      surf(reshape(V(:,s,:),nx,nz)');
   elseif strcmp(d,'x')
      surf(reshape(V(s,:,:),ny,nz)');
   end
end

end
