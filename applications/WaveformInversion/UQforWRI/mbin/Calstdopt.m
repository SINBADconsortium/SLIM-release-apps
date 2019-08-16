expdir = '/scratch/zfang/Result/UQWRI/BGModel/Journal';

n               = [205,451];
A(prod(n),1000) = 0;
k               = 0;

for i = 1:1000
    filename = sprintf([expdir '/NoiseTest%0.4d/m_10.mat'], i);
    if exist(filename,'file');
        k = k +1;
        A(:,k) = vec(ReadAllData(filename));
    end
end
A = A(:,1:k);

for i = 1:prod(n);
    stdv(i)  = std(A(i,:));
    Meanv(i) = mean(A(i,:));
end
stdv  = reshape(stdv,n);
Meanv = reshape(Meanv,n);

figure;imagesc(stdv);colormap redblue
figure;imagesc(Meanv);colormap jet;caxis([1.5,4.5])
