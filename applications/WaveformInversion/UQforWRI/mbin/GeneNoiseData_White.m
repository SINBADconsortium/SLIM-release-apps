[A n d o] = ReadAllData('/scratch/zfang/Data/BG2D/BGDataForSEG2016_15Hz_02s_tristan.mat');
B = squeeze(A);
sigma = .5;
nB        = size(B);
nsub    = nB(1:2);
C     = B;

for i = 1:size(B,3)
        
        rns      = randn(nsub) + sqrt(-1)*randn(nsub);
        beta     = norm(vec(B(:,:,1))) / norm(rns(:));
        C(:,:,i) = B(:,:,i) + rns * beta * sigma;
        Sigma(i) = beta * sigma;
        
end
%    WriteAllData(['./DataForUQPaperBG_White/BGDataNoise_35' num2str(j) '.mat'],reshape(C,n),n,d,o);
WriteAllData(['/scratch/zfang/Data/BG2D/DataSigma_velprior.mat'],Sigma,[size(B,3),1],[d(5),1],[o(5),1]);
