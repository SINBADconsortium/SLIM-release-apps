function show_model(data, colormap_choice, dz, dx)
% syntax: show_model(data, colormap_choice, dz, dx)

if nargin < 2
    colormap_choice = 'gray';
    dz = 1;
    dx = 1;
    physical_distance = 0;
elseif nargin <3
    dz = 1;
    dx = 1;
    physical_distance = 0;
elseif nargin < 4
    dx = dz;
    physical_distance = 1;
elseif nargin == 4
    physical_distance = 1;
end

if ~isreal(data)
    disp('Model contain complex numbers. Only real part is shown')
    data = real(data);
end

if size(data,2)~=1
    disp('Showing a velocity or density model...')
    if physical_distance == 0
        x_axis = (1:size(data,2));
        z_axis = (1:size(data,1));
    else
        x_axis = ((1:size(data,2))-1)*dx;
        z_axis = ((1:size(data,1))-1)*dz;
    end
    figure;
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
    if physical_distance == 0
        xlabel('Horiztontal coordinates');
        ylabel('Depth coordinates');
    else
        xlabel('Horizontal distance (m)')
        ylabel('Depth (m)')
    end
    title('Model')
    pbaspect([dx*size(data,2) dz*size(data,1) 1]);
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal')
    set(get(gca,'XLabel'),'FontName','Arial','FontSize',12,'FontWeight','normal')
    set(get(gca,'YLabel'),'FontName','Arial','FontSize',12,'FontWeight','normal')
    set(get(gca,'Title'),'FontName','Arial','FontSize',12,'FontWeight','normal')
else
    error('Model must be a 2D array.')
end