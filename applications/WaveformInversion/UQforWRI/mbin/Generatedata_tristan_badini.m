%addpath('/Volumes/Users/zfang/Documents/toolsfang/SoftRelease/tools/algorithms/2DFreqModeling');
%IniVel   = '/Volumes/Users/zfang/Documents/ZhilongFang/zfang/Model/BGModel/TrueModel/BG450_tristan.rsf';
addpath(genpath(pwd));
IniVel   = '/scratch/slim/zfang/Model/BG2D/BGIniLinear17.mat';
filename = 'BGDataForPaper2017badIni.mat';

expdir = '/scratch/slim/zfang/Data/BG2D/UQ2017';
curdir = pwd;

if ~exist(expdir,'dir')
        mkdir(expdir)
end
Filename = [expdir '/' filename];
Cfile   = mfilename('fullpath');
Cfile   = [Cfile '.m'];
copyfile(Cfile,expdir);

%expdir = '/Volumes/Users/zfang/Documents/ZhilongFang/zfang/Model/Layer/Data3_reflect';


%if ~exist(expdir,'dir')
%    mkdir(expdir);
%end
%copyfile('Generatedata_tristan.m',expdir)
%cd(expdir)

[v,n,d,o] = ReadAllData(IniVel);

model.n    = n;
model.d    = d;
model.o    = o;
model.nb   = [20,20];
model.freq = [7:9];
model.f0   = 20;
model.t0   = 0.5;
model.xsrc = [0:50:4500];
model.zsrc = [20];
model.xrec = [0:10:4500];
model.zrec = 20;
model.fwave = ones(length(model.freq),1)';
Q          = speye(length(model.xsrc)*length(model.zsrc));

m          = v(:);
m          = 1./m.^2;
Utr        = F(m,Q,model);
Utr        = gather(Utr);

o          = [model.zrec(1),model.xrec(1),model.zsrc(1),model.xsrc(1),model.freq(1)];
if length(model.zrec) == 1
    d(1) = 0;
else
    d(1) = model.zrec(2)-model.zrec(1);
end

if length(model.xrec) == 1
    d(2) = 0;
else
    d(2) = model.xrec(2)-model.xrec(1);
end

if length(model.zsrc) == 1
    d(3) = 0;
else
    d(3) = model.zsrc(2)-model.zsrc(1);
end

if length(model.xsrc) == 1
    d(4) = 0;
else
    d(4) = model.xsrc(2)-model.xsrc(1);
end

if length(model.freq) == 1
    d(5) = 0;
else
    d(5) = model.freq(2)-model.freq(1);
end
label = {'zrec','xrec','zsrc','xsrc','freq'};
ndata(1) = length(model.zrec);
ndata(2) = length(model.xrec);
ndata(3) = length(model.zsrc);
ndata(4) = length(model.xsrc);
ndata(5) = length(model.freq);

WriteAllData(Filename,reshape(Utr,ndata),ndata,d,o);
