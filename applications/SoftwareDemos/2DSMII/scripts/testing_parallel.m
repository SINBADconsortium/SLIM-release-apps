% This script computes data for a 1204x404 model for 10 sources and 64
% frequencies and writes the total CPU time to a file. 
% Run with different number of workers to obtain a scalability table.

setpaths;
curdir = pwd;
expdir = [resultsdir '/testing'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end

% read model
[v,o,d,n]  = odnread([datadir '/marm_vp.odn']);

% model grid
model.o = o;
model.d = d;
model.n = n;

% absorbing boundary
model.nb = [50 50];

% source/receiver grid
model.zsrc = 25;
model.xsrc = linspace(0,11000,10);
model.zrec = 15;
model.xrec = 0:100:11000;

% frequencies
model.freq = 4.25:.25:20;

% wavelet
model.f0 = 10;
model.t0 = 0;

% source matrix
Q = speye(length(model.xsrc));

spmd,
    np = numlabs;
end
np = np{1};

for k = 1:5
    tic;
    D = F(1e6./v.^2,Q,model);
    T(k)=toc;    
end


fid = fopen([resultsdir '/testing/time.dat'],'a');

fclose(fid);
