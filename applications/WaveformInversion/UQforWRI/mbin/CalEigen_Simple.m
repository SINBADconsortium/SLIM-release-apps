vfile     = '/scratch/zfang/Model/SimpleLayer/SimpleLayerkm.mat';
addpath(genpath(pwd));

nb        = 20;

[vel n d o]     = ReadAllData(vfile);
vel             = vel;
m               = 1./vel(:).^2;
model.o  = o;
model.d  = d;
model.n  = n;
model.nb = nb*[1 1];
model.freq = 2;


model.xrec = [0:50:3000];
model.zrec = 50;

[lambda, A ,Pr] = CalculateEigenvalue(model,m);
