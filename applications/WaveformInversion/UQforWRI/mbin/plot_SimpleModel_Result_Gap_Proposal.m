ResDir = '/scratch/zfang/Result/UQWRI/SimpleLayer/Journal/RTO_SMP_Cent6';
ResDir2 = '/scratch/zfang/Result/UQWRI/SimpleLayer/Journal/RTO_SMP_GAP';
FigDir = '/scratch/zfang/Figure/UQWRI/Paper/Journal/ForProposal';
addpath('/home/zfang/Documents/Project/zfang-tool');
addpath(genpath(pwd));
Cfile   = mfilename('fullpath');
Cfile   = [Cfile '.m'];
copyfile(Cfile,FigDir);

%% File name

STD_SMP_RTO_FILE  = [ResDir '/STD_RTO.mat'];
MEAN_SMP_RTO_FILE = [ResDir '/MEAN_RTO.mat'];
SMP_RTO_FILE      = [ResDir '/SMP_RTO.mat'];
STD_SMP_INV_FILE  = [ResDir2 '/STD_RTO.mat'];
MEAN_SMP_INV_FILE = [ResDir2 '/MEAN_RTO.mat'];
SMP_INV_FILE      = [ResDir2 '/SMP_RTO.mat'];
STD_CHOL_FILE     = [ResDir '/STD_CHOL.mat'];
MAP_FILE          = ['/scratch/zfang/Model/SimpleLayer/SimpleLayerkm.mat'];
MTRUE_FILE        = ['/scratch/zfang/Model/SimpleLayer/SimpleLayerkm.mat'];
MINI_FILE         = ['/scratch/zfang/Model/SimpleLayer/Ini_Vel.mat'];

%% Read file

[m_MAP, n, d, o]   = ReadAllData(MAP_FILE);
m_true             = ReadAllData(MTRUE_FILE);
m_ini              = ReadAllData(MINI_FILE);
[STD_RTO]          = ReadAllData(STD_SMP_RTO_FILE);
[STD_INV]          = ReadAllData(STD_SMP_INV_FILE);
[STD_CHOL]         = ReadAllData(STD_CHOL_FILE);
[MEAN_RTO]         = ReadAllData(MEAN_SMP_RTO_FILE);
[MEAN_INV]         = ReadAllData(MEAN_SMP_INV_FILE);
% [SMP_RTO]          = ReadAllData(SMP_RTO_FILE);
% [SMP_INV]          = ReadAllData(SMP_INV_FILE);

%% Calculate std of prior
nlayercst       = 5;
SIGMA_mp_ratior = 0.1;
velp            = m_true;
SIGMA_z         = [ones(nlayercst,1)*1e-1; vec(velp(nlayercst+1:n(1))) * SIGMA_mp_ratior];
SIGMA_mp        = repmat(SIGMA_z,1,n(2));

%% Flag of figure 

FlagModel = 1;
FlagCross = 1;
FlagCross2 = 0;


% parameters for plotting Model
str_Model   = {'Layer_True','Layer_Ini','Layer_Inistd','Layer_MAP','Layer_RTOMEAN','Layer_RTOSTD','Layer_GAPMEAN','Layer_GAPSTD','Layer_CHOLSTD'};
M_Model     = {m_true,m_ini,SIGMA_mp,m_MAP, MEAN_RTO, STD_RTO, MEAN_INV, STD_INV, STD_CHOL};
C_MAP       = {jet,   jet,  redblue, jet,   jet,      redblue, jet,      redblue, redblue};
cax_M       = [1.5 2.5];
cax_STD     = [0.003 0.25];
CAX         = {cax_M, cax_M,cax_STD, cax_M, cax_M,     cax_STD, cax_M,   cax_STD, cax_STD};
Paper_Model = [10 4];

% parameters for plotting cross line
xpos           = [ 500 1500 2500];
Paper_STDcmp   = [5,8];

if FlagModel > 0
    [zz xx] = odn2grid(o,d,n);
    for i = 1:length(str_Model)
        filename = [FigDir '/' str_Model{i}];
        mtmp     = M_Model{i};
        CMAPtmp  = C_MAP{i};
        mtmp(1,:)   = mtmp(2,:);
        mtmp(end,:) = mtmp(end-1,:);
        mtmp(:,1)   = mtmp(:,2);
        mtmp(:,end) = mtmp(:,end-1);
        figure;fidx=gcf;
        imagesc(xx,zz,mtmp);colormap(CMAPtmp);
        xlabel('Lateral [m]','fontsize',20);ylabel('Depth [m]','fontsize',20);
        hc = colorbar; title(hc,'km/s');caxis(CAX{i});
        set(gca,'fontsize',20);
        PrintFigure(fidx,Paper_Model,filename);pause(0.5)
        convert_command = ['convert ' filename '.pdf ' filename '.png'];
        status          = unix(convert_command);
    end
end

if FlagCross > 0
    [zz xx] = odn2grid(o,d,n);
    I   = floor((xpos-o(2))/d(2)) + 1;
    for i = 1:length(I)
        velRTO = MEAN_RTO(:,I(i));
        velINV = MEAN_INV(:,I(i));
        stdRTO = STD_RTO(:,I(i));
        stdINV = STD_INV(:,I(i));
        stdRTO(1)   = stdRTO(2);
        stdRTO(end) = stdRTO(end-1);
        stdINV(1)   = stdINV(2);
        stdINV(end) = stdINV(end-1);
        figure;fidx=gcf;
        plot_MeanCmp(zz,velRTO,stdRTO,velINV,stdINV);
        xlabel('Velocity [km/s]','fontsize',20)
        ylabel('Depth [m]','fontsize',20)
        h   = legend('Without Gap','With Gap');
        set(h,'fontsize',14);set(gca,'fontsize',20);
        filename = [FigDir '/Layer_STDcmpGAPX_' num2str(xpos(i))];
        PrintFigure(fidx,Paper_STDcmp,filename);pause(0.5)
        convert_command = ['convert ' filename '.pdf ' filename '.png'];
        status          = unix(convert_command);
    end
end

if FlagCross2 > 0
    [zz xx] = odn2grid(o,d,n);
    I   = floor((xpos-o(2))/d(2)) + 1;
    for i = 1:length(I)
        velRTO = MEAN_RTO(:,I(i));
        velINV = MEAN_RTO(:,I(i));
        stdRTO = STD_RTO(:,I(i));
        stdINV = SIGMA_mp(:,I(i));
        stdRTO(1)   = stdRTO(2);
        stdRTO(end) = stdRTO(end-1);
        stdINV(1)   = stdINV(2);
        stdINV(end) = stdINV(end-1);
        figure;fidx=gcf;
        plot_MeanCmp(zz,velRTO,stdRTO,velINV,stdINV);
        xlabel('Velocity [km/s]','fontsize',20)
        ylabel('Depth [m]','fontsize',20)
        h   = legend('Inversion result','Prior model');
        set(h,'fontsize',14);set(gca,'fontsize',20);
        filename = [FigDir '/Layer_STDcmpPriorX_' num2str(xpos(i))];
        PrintFigure(fidx,Paper_STDcmp,filename);pause(0.5)
        convert_command = ['convert ' filename '.pdf ' filename '.png'];
        status          = unix(convert_command);
    end
end


