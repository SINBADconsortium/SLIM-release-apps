function show_model(data, options)
% Syntax:
% show_model(data, options)
%
% Note:
% This function is called by Migration_with_SRM. Not designed as a stand-alone
% function.
% 
% Description:
% Display 2D velocity/density models
%
% Input list:
% 1. data: apprarently is the model you want to show :D
% 2. options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
%
% Output list: None. Figures will be displayed on screen.
%
% By Ning Tu @ SLIM, UBC
% Date: Dec/14/2012

disp('Showing model...You need to make sure that the model is of correct dimension.')

colormap_choice = options.color_model;
if isempty(colormap_choice)
    colormap_choice = 'gray';
end

if ~isreal(data)
    disp('Model data contain complex numbers. Only real part is shown')
    data = real(data);
end

if size(data,2)~=1
    disp('Showing a velocity or density model...')
    x_axis = ([1:size(data,2)]-1)*options.dx;
    z_axis = ([1:size(data,1)]-1)*options.dz;
    imagesc(x_axis,z_axis,data)
    colormap(colormap_choice)
    val_min = min(data(:));
    val_max = max(data(:));
    val_limit = max(abs(val_min),abs(val_max));
    if val_min < 0
        caxis([-0.95*val_limit, 0.95*val_limit])
    else
        caxis([0.99*val_min, 1.01*val_max]);
    end
    pbaspect([options.dx*size(data,2) options.dz*size(data,1) 1]);
    xlabel(['Horiztontal distance (m)']);
    ylabel(['Depth (m)']);
    title('Model')
else
    error('Model must be a 2D array.')
end