function PrintFigure(fidx, Papersize, strfile)
	%% PrintFigure(fidx, Papersize)
	%  fidx  - index of figure
	%  Papersize - paper size
	%  strfile   - file name
	
    figure(fidx);
	set(gcf, 'PaperSize', Papersize);
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperPosition', [0 0 Papersize(1) Papersize(2)]);
	set(gcf, 'renderer', 'painters');
%    set(fidx, 'Renderer', 'opengl');
	print(fidx,'-depsc',[strfile '.eps'],'-loose');
	scommand                = ['unset DYLD_LIBRARY_PATH;convert ' strfile '.eps ' strfile '.png'];
	[status,result]=system(scommand);
        scommand                = ['epstopdf ' strfile '.eps'];
        [status,result]=system(scommand); 
