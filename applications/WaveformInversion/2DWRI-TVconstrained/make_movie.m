%make a movie

figure(99);
clf;
numframes = model.batches+1;
velocity_movie(numframes).cdata = [];
velocity_movie(numframes).colormap = [];
for f = 1:numframes
    
    imagesc(model.xt,model.zt,1./sqrt(mb(:,:,f)));
    title(['velocity estimate after ' num2str(f-1) ' frequency batches']);
    colorbar; caxis([min(model.vmin(:)) max(model.vmax(:))]);
    
    velocity_movie(f) = getframe(gcf);
    
end

%movie(99,velocity_movie,1,4)