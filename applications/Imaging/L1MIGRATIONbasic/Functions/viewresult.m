function viewresult(name,n,m)


if nargin > 0 
    result = load('./result/',name);
    if nargin > 1
        result = reshape(result,n,m);
    end
    imagesc(result);
end

% mamousi model FWI result with wavelet


load ./result/Mmodel_result_wavtest10
figure;imagesc(vel1)
caxis([1500 5500]);colorbar
title('Mamousi model FWI result with source wavelet without stack')

load ./result/BPmodel_P4result10 
figure;imagesc(vel1)
caxis([1500 5500]);colorbar
title('Mamousi model FWI result with source wavelet with stack')

% BP FWI result without wavelet

load ./result/BPmodel_F_result10 
figure;imagesc(vel1)
caxis([1500 5500]);colorbar
title('BP model FWI result without source wavelet')

% Mamousi model RTM with L2 same downsampling

load ./result/Mmodel_RTML2asL1_update1
figure;imagesc(reshape(dm,128,384))
title('Mamousi model Linearized migration with L2')


% Mamousi model RTM with L1 without renewal


load ./result/Mmodel_RTMworn_update1
figure;imagesc(reshape(dm(:,end),128,384))
title('Mamousi model Linearized migration with L1 without renewals')

% Mamousi model RTM with L1 with renewal


load ./result/Mmodel_RTMwrn_update_RSRF1
figure;imagesc(reshape(dm(:,end),128,384))
title('Mamousi model Linearized migration with L1 with renewals')