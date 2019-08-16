ResDir = '/scratch/zfang/Result/UQWRI/BGModel/Journal_1e3_02/NoiseTest/RTO_SMP4';
OutputDir = '/scratch/zfang/Result/UQWRI/BGModel/Journal_1e3_02/NoiseTest';

alpha  = 0.95;
n = [205,451];
V = [];
for i = 1:200
    str = [ResDir '/Smp_' num2str(i) '.mat'];
    if exist(str,'file')
        load(str);
        V =[V Smp];
    end
end

for i = 1:size(V,1)
    stdv(i)  = std(V(i,:));
    CI(i,:)  = confidence_interval(V(i,:),alpha);
end

CIout(:,:,1) = reshape(CI(:,1),n);
CIout(:,:,2) = reshape(CI(:,2),n);

outputfile = [OutputDir '/CI_RTO.mat'];
WriteAllData(outputfile, CIout, [n(1) n(2) 2], [10 10 1], [0 0 0]);
