ResDir1 = '/scratch/zfang/Result/UQWRI/BGModel/Journal_1e3_02/NoiseTest';
ResDir2 = '/scratch/zfang/Result/UQWRI/BGModel/Journal_Test_PriorIni_1000_2';
FigDir  = '/scratch/zfang/Figure/UQWRI/Paper/Journal/SEG';
addpath(genpath(pwd));

xpos           = [1000 2500 4000];
Paper_STDcmp   = [5,8];

xid = [101 251 401];
sid = [10:30:100];

MEANfile = [ResDir1 '/m_10.mat'];
STDfile  = [ResDir1 '/STD1e3.mat'];
CIfile   = [ResDir1 '/CI_RTO.mat'];

[m n d o]    = ReadAllData(MEANfile);
stdv = ReadAllData(STDfile);
CI   = ReadAllData(CIfile);

[z x] = odn2grid(o,d,n);

for i = 1:length(xid)
    k    = xid(i);
    mtmp = m(:,k);
    stdtmp = stdv(:,k);
    CI1    = mtmp - CI(:,k,1);
    CI2    = -mtmp + CI(:,k,2);
    
    for j=1:length(sid)
        str = sprintf('/NoiseTest0%0.3dTest_Test10/m_10.mat',sid(j));
        str = [ResDir2 str];
        if exist(str,'file')
            A = ReadAllData(str);
%            V = [V real(A(:,k))];
        end
        figure;fidx=gcf;
        [l p]=boundedline(z,mtmp,[CI1 CI2]);
        ylabel('Velocity [km/s]','fontsize',20)
        xlabel('Depth [m]','fontsize',20)
        h   = legend([p(1)],'0.95 Confidence interval');
         hold on
        V = [];
    
    
   
%        MEANV(:,i) = mean(V,2);
%        for l = 1:n(1)
%                STDV(l,i) =std(V(l,:));
%        end
        plot(z,A(:,k),'red','linewidth',2); view(90,90);ylim([0,6]);xlim([0,2050]);
    
        set(h,'fontsize',14);set(gca,'fontsize',20);
        filename = [FigDir '/BG_CIX_' num2str(xpos(i)) '_' num2str(sid(j))];
        PrintFigure(fidx,Paper_STDcmp,filename);pause(0.5)
        convert_command = ['convert ' filename '.pdf ' filename '.png'];
        status          = unix(convert_command);
    end
%    figure;plot_MeanCmp(z,mtmp,stdtmp,MEANV(:,i),STDV(:,i));view(90,90);
end

