function plot_vel_well(data1,data2,weight,x,t,xp,tp,CAXIS,CMAP)
	%  plot_wiggle(data1,data2,nskip,weight,x,t)
	%
	%
	% data1: background data
	% data2: wiggle data, same size as data1
	% nskip: plot wiggle every 'nskip' traces
	% weight: weight data2, if the amp is very small;
	% x,t: x axis and t axis;
	%
	% %

	
    imagesc(x,t,data1);caxis(CAXIS);colormap(CMAP)
    
    idx = find(x>xp,1);
    m   = idx -1;
    
        
    line(weight.*data2+x(m),tp,'color','k');
        	
