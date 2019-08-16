function PrintWebFigure(v,n,d,o,cmap,Caxis,filename,Pagesize)
        if nargin < 8
                Pagesize = [12 4];
        end
	curdir = pwd;
	cd('/home/slim/zfang/public_html/Figure/Chevron/tmp');
	plotModel(v,n,d/1000,o);set(gca,'fontsize',20);
	xlabel('Lateral [km]');ylabel('Depth [km]');
	colormap(cmap);colorbar;caxis(Caxis);
	fidx = gcf;
	PrintFigure(gcf,Pagesize,filename);
	!rm *.eps
        !rm *.pdf
	cd(curdir)
