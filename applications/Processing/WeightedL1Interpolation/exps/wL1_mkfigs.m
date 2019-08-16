% Plotting interpolation results using weighted L1 minimization of a
% subsampled seismic line from the Gulf of Suez.
% Source/receiver spacing is 12.5m
% Temporal resolution is 4ms.
% Iime variable selects an arbitrary time slice from the seismic line to display the subsampling and results
% So it can take any value between 1 and dim(1), where dim = size(D);

% Copyright 2013 Hassan Mansour (hassanm@cs.ubc.ca)


close all;


% choose data slice 
shot_num = 84;
time = 350;

%% Source receiver domain
% set directories
label = 'wL1FreqSR';
outputdir = ['../results/' label];
if ~exist(outputdir,'dir')
    error('wL1_freq_SR.m, wL1_freq_MH.m, and wL1_offset_TM.m must be run before plotting figures.');
end
load([outputdir, '/wL1_freq_SR.mat']);

dim = size(D);
DSG = reshape(D(:,:,shot_num),dim(1),dim(2));
Dt = reshape(D(time,:,:),dim(2), dim(3));

% Original data and subsampling
figure; imagesc(Dt); colormap(seiscol); caxis([-50, 50]); 
title('Original time slice'); xlabel('Source number'); ylabel('Receiver number');
figure; imagesc(mask); colormap(gray); 
title('Subsampling mask'); xlabel('Source number'); ylabel('Receiver number');
figure; imagesc(Dt.*mask); colormap(seiscol); caxis([-50, 50])
title('Subsampled time slice'); xlabel('Source number'); ylabel('Receiver number');

% Fully sampled and reconstructed shot gather
Dest = real(Dest);

DestSG = reshape(Dest(:,:,shot_num),dim(1),dim(2));
%clear Dest

figure;imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3, DSG);colormap(seiscol); caxis([-50 50])
title('Original shot gather'); xlabel('Distance (m)'); ylabel('Time (sec)');
figure;imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3, (DestSG));colormap(seiscol); caxis([-50 50])
title('Weighted L_1 minimization in SR'); xlabel('Distance (m)'); ylabel('Time (sec)');

Err_wl1SR = DSG-DestSG;
figure;imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,(Err_wl1SR));colormap(seiscol); caxis([-50, 50]) 
title('Weighted L_1 error image'); xlabel('Distance (m)'); ylabel('Time (sec)');

% compute shot gather relative error
for s = 1:dim(2)
    EwL1freqSR_source(s) = norm(vec(D(:,:,s)) - vec(Dest(:,:,s)))/norm(vec(D(:,:,s)));
end


%% Midpoint offset domain
label = 'wL1FreqMH';
outputdir = ['../results/' label];
if ~exist(outputdir,'dir')
    error('wL1_freq_SR.m, wL1_freq_MH.m, and wL1_offset_TM.m must be run before plotting figures.');
end
load([outputdir, '/wL1_freq_MH.mat']);

dim = size(DH);
DHt = reshape(DH(time,:,:),dim(2), dim(3));
clear DH

% Original data and subsampling
figure; imagesc([-dim(2):dim(2)], [1:dim(2)]*12.5,DHt); colormap(seiscol); caxis([-50, 50])
title('Original time slice'); xlabel('Offset number'); ylabel('Distance (m)');
figure; imagesc([-dim(2):dim(2)], [1:dim(2)]*12.5,maskH); colormap(gray); 
title('Subsampling mask'); xlabel('Offset number'); ylabel('Distance (m)');
figure; imagesc([-dim(2):dim(2)], [1:dim(2)]*12.5,DHt.*maskH); colormap(seiscol); caxis([-50, 50])
title('Subsampled time slice'); xlabel('Offset number'); ylabel('Distance (m)');

% Fully sampled and reconstructed shot gather
dim = size(Dest);
DestSG = reshape(Dest(:,:,shot_num),dim(1),dim(2));
%clear Dest

figure;
imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,DestSG);colormap(seiscol); caxis([-50 50])
title('Weighted L_1 minimization in MH'); xlabel('Distance (m)'); ylabel('Time (sec)');

Err_wl1MH = DSG - DestSG;
figure; imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,(Err_wl1MH));colormap(seiscol); caxis([-50 50]) 
title('Weighted L_1 in MH error image'); xlabel('Distance (m)'); ylabel('Time (sec)');

% compute shot gather relative error
for s = 1:dim(2)
    EwL1freqMH_source(s) = norm(vec(D(:,:,s)) - vec(Dest(:,:,s)))/norm(vec(D(:,:,s)));
end



%% Time midpoint domain
label = 'wL1OffsetTM';
outputdir = ['../results/' label];
if ~exist(outputdir,'dir')
    error('wL1_freq_SR.m, wL1_freq_MH.m, and wL1_offset_TM.m must be run before plotting figures.');
end
load([outputdir, '/wL1_offset_TM.mat']);

maskH_vol = permute(repmat(maskH, [1,1,dim(1)]),[3 1 2]);


% Original, subsampled, and reconstructed near offset
offset = dim(2)+5;
DHoff = reshape(DH(:,:,offset),dim(1), dim(2));
figure; imagesc( [1:dim(2)]*12.5,[1:dim(1)]*4e-3, DHoff);
colormap(seiscol); caxis([-50, 50]);
title('Original +5-offset gather'); xlabel('Distance (m)'); ylabel('Time (sec)');

figure; imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,DHoff.*maskH_vol(:,:,offset));
colormap(seiscol); caxis([-50 50])
title('Subsampled shot gather'); xlabel('Distance (m)'); ylabel('Time (sec)');

DHwl1off = reshape(real(DHest(:,:,offset)),dim(1), dim(2));
figure; imagesc( [1:dim(2)]*12.5,[1:dim(1)]*4e-3, DHwl1off);
colormap(seiscol); caxis([-50, 50]);
title('Weighted L_1 recovered +5-offset'); xlabel('Distance (m)'); ylabel('Time (sec)');

ErrHwl1off = DHoff - DHwl1off;
figure; imagesc( [1:dim(2)]*12.5,[1:dim(1)]*4e-3, ErrHwl1off);title('Weighted L_1 +5-offset error image'); 
colormap(seiscol); caxis([-50, 50]);
title('Weighted L_1 +5-offset error image'); xlabel('Distance (m)'); ylabel('Time (sec)');

% Fully sampled and reconstructed shot gather
DestSG = reshape(Dest(:,:,shot_num),dim(1),dim(2));
%clear Dest

%figure; imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,bSG);colormap(seiscol); caxis([-50 50])
title('Subsampled shot gather'); xlabel('Distance (m)'); ylabel('Time (sec)');

figure;
imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,DestSG);colormap(seiscol); caxis([-50 50])
title('Weighted L_1 in Offset'); xlabel('Distance (m)'); ylabel('Time (sec)');

Err_wl1Offset = DSG - DestSG;

figure;imagesc([1:dim(2)]*12.5, [1:dim(1)]*4e-3,(Err_wl1Offset));colormap(seiscol); caxis([-50, 50]) 
title('Weighted L_1 in Offset error image'); xlabel('Distance (m)'); ylabel('Time (sec)');

% compute shot gather relative error
for s = 1:dim(2)
    EwL1OffsetTM_source(s) = norm(vec(D(:,:,s)) - vec(Dest(:,:,s)))/norm(vec(D(:,:,s)));
end

%% compare the reconstruction SNRs of all shot gathers for the three interpolation schemes
figure; plot(1:dim(2), -20*log10(EwL1freqSR_source), 1:dim(2), -20*log10(EwL1freqMH_source),...
    1:dim(2), -20*log10(EwL1OffsetTM_source)); 
legend('W-L_1 on freq. in SR', 'W-L_1 on freq. in MH','W-L_1 on offsets in TM'); 
xlabel('Shot gather number'); ylabel('SNR')
