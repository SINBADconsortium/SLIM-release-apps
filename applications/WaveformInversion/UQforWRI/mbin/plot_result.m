ResDir = '/scratch/zfang/Result/UQWRI/BGModel/Journal/FinalResult';
FigDir = '/scratch/zfang/Figure/UQWRI/Paper/Journal';
addpath('/home/zfang/Documents/Project/zfang-tool');
addpath(genpath(pwd));
Cfile   = mfilename('fullpath');
Cfile   = [Cfile '.m'];
copyfile(Cfile,FigDir);

%% File name

STD_SMP_RTO_FILE  = [ResDir '/STD_SMP_RTO.mat'];
MEAN_SMP_RTO_FILE = [ResDir '/MEAN_SMP_RTO.mat'];
SMP_RTO_FILE      = [ResDir '/SMP_RTO.mat'];
STD_SMP_INV_FILE  = [ResDir '/STD_SMP_INV.mat'];
MEAN_SMP_INV_FILE = [ResDir '/MEAN_SMP_INV.mat'];
SMP_INV_FILE      = [ResDir '/SMP_INV.mat'];
MAP_FILE          = ['/scratch/zfang/Result/UQWRI/BGModel/Journal/NoiseTest/m_10.mat'];
MTRUE_FILE        = ['/scratch/zfang/Model/BG2D/BG450kmTrue.mat'];
MINI_FILE         = ['/scratch/zfang/Model/BG2D/Ini_Vel.mat'];

%% Read file

[m_MAP, n, d, o]   = ReadAllData(MAP_FILE);
m_true             = ReadAllData(MTRUE_FILE);
m_ini              = ReadAllData(MINI_FILE);
[STD_RTO]          = ReadAllData(STD_SMP_RTO_FILE);
[STD_INV]          = ReadAllData(STD_SMP_INV_FILE);
[MEAN_RTO]         = ReadAllData(MEAN_SMP_RTO_FILE);
[MEAN_INV]         = ReadAllData(MEAN_SMP_INV_FILE);
[SMP_RTO]          = ReadAllData(SMP_RTO_FILE);
[SMP_INV]          = ReadAllData(SMP_INV_FILE);

%% Calculate std of prior
nlayercst       = 5;
SIGMA_mp_ratior = 0.3;
velp            = m_ini;
SIGMA_z         = [ones(nlayercst,1)*1e-1; vec(velp(nlayercst+1:n(1))) * SIGMA_mp_ratior];
SIGMA_mp        = repmat(SIGMA_z,1,n(2));

%% Flag of figure 

FlagModel = 0;
FlagCross = 1;


% parameters for plotting Model
str_Model   = {'True','Ini','Inistd','MAP','RTOMEAN','RTOSTD','INVMEAN','INVSTD'};
M_Model     = {m_true,m_ini,SIGMA_mp,m_MAP, MEAN_RTO, STD_RTO, MEAN_INV, STD_INV};
C_MAP       = {jet,   jet,  redblue, jet,   jet,      redblue, jet,      redblue};
cax_M       = [1.5 4.5];
cax_STD     = [0.0074 0.1545];
CAX         = {cax_M, cax_M,[0.0074,1.5], cax_M, cax_M,     cax_STD, cax_M,   cax_STD};
Paper_Model = [10 4];

% parameters for plotting cross line
xpos           = [ 1000 2500 4000];
Paper_STDcmp   = [5,8];

if FlagModel > 0
    [zz xx] = odn2grid(o,d,n);
    for i = 1:length(str_Model)
        filename = [FigDir '/BG' str_Model{i}];
        mtmp     = M_Model{i};
        CMAPtmp  = C_MAP{i};
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
    I   = floor((xpos-o(2))/d(2)) + 1;
    for i = 1:length(I)
        velRTO = MEAN_RTO(:,I(i));
        velINV = MEAN_INV(:,I(i));
        stdRTO = STD_RTO(:,I(i));
        stdINV = STD_INV(:,I(i));
        figure;fidx=gcf;
        plot_MeanCmp(zz,velRTO,stdRTO,velINV,stdINV);
        xlabel('Velocity [km/s]','fontsize',20)
        ylabel('Depth [m]','fontsize',20)
        h   = legend('RTO','Inversion');
        set(h,'fontsize',14);set(gca,'fontsize',20);
        filename = [FigDir '/BG_STDcmpX_' num2str(xpos(i))];
        PrintFigure(fidx,Paper_STDcmp,filename);pause(0.5)
        convert_command = ['convert ' filename '.pdf ' filename '.png'];
        status          = unix(convert_command);
    end
end










