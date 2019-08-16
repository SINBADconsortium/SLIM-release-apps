ResDir1 = '/scratch/zfang/Result/UQWRI/BGModel/Journal/FinalResult_20';
ResDir2 = '/scratch/zfang/Result/UQWRI/BGModel/Journal_Test_PriorIni_1000_2';

sid = [1:1:1000];

V = [];
for j=1:length(sid)
    str = sprintf('/NoiseTest0%0.3dTest_Test10/m_10.mat',sid(j));
    str = [ResDir2 str];
    if exist(str,'file')
        A = ReadAllData(str);
        V = [V real(A(:))];
    end
    
end
MEANV = mean(V,2);
for l = 1:size(V,1)
    STDV(l) =std(V(l,:));
end

MEANV = reshape(MEANV,205,451);
STDV = reshape(STDV,205,451);

str1 = [ResDir1 '/STD_SMP_INV.mat'];
str2 = [ResDir1 '/MEAN_SMP_INV.mat'];

WriteAllData(str1,STDV,[205,451],[10 10],[0 0]);
WriteAllData(str2,MEANV,[205,451],[10 10],[0 0]);