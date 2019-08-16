function MakeMovie_Long(Data, filename, DataSize, FigureSize, fcaxis, FrameRate, Cmap)
%% Function to make a movie with .avi format
% 	MakeMovie_Long(Data, filename, DataSize, PaperSize, fcaxis, FrameRate)
%   Input: 
%   Data      : Data = [d1 d2 d3 ... dn];
%   filename  : Output filename should be ended with '.avi'
%   DataSize  : The true 2D size of your data
%   PaperSize : The size of your movie, PaperSize = [Width, Height]
%   fcaxis    : color axis for your movie
%   FrameRate : 
	

model.n   = DataSize;
aviname   = filename;
writerobj = VideoWriter(aviname);
writerobj.FrameRate = FrameRate; 
Width     = FigureSize(1);
Height    = FigureSize(2);
open(writerobj)
for i = 1:size(Data,2)
	hFig = figure;
	set(hFig, 'Position', [600 300 Width Height])
    imagesc(reshape(Data(:,i),model.n));caxis(fcaxis);colormap(Cmap)
    F = getframe;
    movie2(:,:,:,i) = F.cdata;
    writeVideo(writerobj,F.cdata);
    close
    i
end
close(writerobj)
