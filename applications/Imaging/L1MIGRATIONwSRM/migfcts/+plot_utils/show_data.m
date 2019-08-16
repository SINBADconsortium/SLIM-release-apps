function show_data(data,options,compare_traces,offset)
% Syntax:
% show_data(data,options,compare_traces,offset)
%
% Note:
% This function is called by Migration_with_SRM. Not designed as a stand-alone
% function.
%
% Description:
% Display seismic data:
% 1. if data is 3D, will show each 2D slice, with a GUI to scroll over slices
% 2. if data is 2D and "compare_traces" is 0, will display data and axises
% 3. if data is 2D and "compare_traces" is 1, will display each columns
% 4. if data is 1D, will display the trace with a correct axis
%
% Input list:
% 1. data: apprarently is the data you want to show :D
% 2. options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
% 3. compare_traces: optional input. If is 1, will display each column of a 2D
%		data instead of display the 2D data.
% 4. offset: only available for 2D and 1D data display, is a 1x2 vector.
%		offset(1) is the offset in physical distance axis, and offset(2)
%		is the offset in time axis.
%
% Output list: None. Figures will be displayed on screen.
%
% By Ning Tu @ SLIM, UBC
% Date: Feb/14/2012 (Happy Valentine's Day!)

if nargin < 3
    compare_traces = 0;
    offset = [0 0];
elseif nargin < 4
    offset = [0,0];
end

disp('Showing data...You need to make sure that data is of correct dimension.')

colormap_choice = options.color_data;
if isempty(colormap_choice)
    colormap_choice = 'gray';
end

if ~isreal(data)
    disp('Data contain complex numbers. Only real part is shown')
    data = real(data);
end

if size(data,3)~=1
    disp('Showing a data cube.')
    scrollData(data)
    colormap(colormap_choice)
    title('3D data, scroll to see different data slices.')
elseif size(data,2)~=1 && ~compare_traces
    disp('Showing a 2D data slice.')
    x_axis = (options.x_rcv_grid-1+offset(1))*options.dx;
    t_axis = ([1:size(data,1)]-1+offset(2))*options.dt;
    imagesc(x_axis,t_axis,data)
    colormap(colormap_choice)
    val_limit = max(abs(data(:)));
    caxis([-0.95*val_limit, 0.95*val_limit])
    xlabel(['Horiztontal distance (m)']);
    ylabel('Time (S)')
    title('2D data slice.')
else
    disp('Plotting a 1D data trace.')
    t_axis = ([1:size(data,1)]-1+offset(2))*options.dt;
    plot(t_axis,data)
    xlabel('Time(s)')
    title('1D data trace.')
end

set(gca,'FontName','Arial','FontSize',12,'FontWeight','normal')
set(get(gca,'XLabel'),'FontName','Arial','FontSize',12,'FontWeight','normal')
set(get(gca,'YLabel'),'FontName','Arial','FontSize',12,'FontWeight','normal')
set(get(gca,'Title'),'FontName','Arial','FontSize',12,'FontWeight','normal')