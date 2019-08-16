function displayDatacube(data)
% blocks until user closes figure

dims = size(data);
nt   = dims(1);
nr   = dims(2);
ns   = dims(3);

crossplane_t = ceil(nt/3);
crossplane_r = ceil(nr/2);
crossplane_s = ceil(ns/2);

fh = figure;
title('Corssing-plane view of input seismic datacube')
colormap(gray)


% Plot same-time slice
subplot(2,2,1)
imagesc(squeeze(data(crossplane_t,:,:)))
title(['Const-time slice ' num2str(crossplane_t) '/' num2str(nt) ' (should be diagonally symmetric)'])

% Plot same-receiver slice
subplot(2,2,2)
imagesc(squeeze(data(:,crossplane_r,:)))
title(['Const-receiver slice, receiver ' num2str(crossplane_r) '/' num2str(nr)])


% Plot same-shot slice
subplot(2,2,3)
imagesc(squeeze(data(:,:,crossplane_s)))
title(['Const-shot slice, shot ' num2str(crossplane_s) '/' num2str(ns)])


% Plot a trace
subplot(2,2,4)
plot(squeeze(data(:,crossplane_r,crossplane_s)))
title(['Single trace, receiver ' num2str(crossplane_r) ', shot ' num2str(crossplane_s)])

disp('Please review the plotted figure of loaded datacube...')
waitfor(fh)
return