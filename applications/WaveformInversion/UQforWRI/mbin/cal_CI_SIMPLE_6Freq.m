ResDir = '/scratch/zfang/Result/UQWRI/SimpleLayer/Journal_Test_6Freq/NoiseTest102';
OutputDir = '/scratch/zfang/Result/UQWRI/SimpleLayer/Journal/RTO_SMP_6Freq';

alpha  = 0.9;
n = [31,61];
V = [];
for i = 1:200
    str = [ResDir '/Smp' num2str(i) '.mat'];
    if exist(str,'file')
        load(str);
        V =[V VV];
    end
end

for i = 1:size(V,1)
    stdv(i)  = std(V(i,:));
    CI(i,:)  = confidence_interval(V(i,:),alpha);
end

CIout(:,:,1) = reshape(CI(:,1),n);
CIout(:,:,2) = reshape(CI(:,2),n);

outputfile = [OutputDir '/CI_INV.mat'];
WriteAllData(outputfile, CIout, [n(1) n(2) 2], [50 50 1], [0 0 0]);
