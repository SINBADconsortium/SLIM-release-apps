vfile     = '/scratch/slim/zfang/Model/BG2D/Ini_Vel.mat';
addpath(genpath(pwd));

nb        = 20;

[vel n d o]     = ReadAllData(vfile);
vel             = vel;
m               = 1./vel(:).^2;
model.o  = o;
model.d  = d;
model.n  = n;
model.nb = nb*[1 1];
model.freq = 10;


model.xrec = [0:10:4500];
model.zrec = 40;

[lambda, A ,Pr] = CalculateEigenvalue(model,m);