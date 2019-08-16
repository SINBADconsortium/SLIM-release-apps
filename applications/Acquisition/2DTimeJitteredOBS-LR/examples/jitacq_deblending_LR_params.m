% This script sets the parameters for time-jittered (blended) marine acquisition, and conventional data
% recovery by rank-minimization. 
%
% Acquisition scenario: time-jittered OBC marine acquisition with one source vessel and two airgun arrays
% (Acquisition scenario is same as proposed by Haneet Wason)

% Set a label for the experiment
label = 'TimeJitAcq_1boat_2array_LR';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% Input data file
fname_input = [datadir '/GulfOfSuez178.su'];  

% Seismic data parameters
% (nt_org, nr_org, ns_org): time, receiver, and source samples for all the data
% (nt, nr, ns): time, receiver, and source samples for a subset of the data
nt_org = 1024;
nr_org = 178;
ns_org = 178;
nt = 1024;
nr = 128; 
ns = 128; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 25.0;     
ds = 25.0;     

% Parameters for time-jittered acquisition
rndfactor = [1000 100];
p = 2;
nboats = 1;
randseed = [8 3];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = 10.0;
delayboat = 0;
fig = 'no';

% Flag to perform dottest for the sampling operator
test = 0;
iteration = 200;
initiweight = 1e-6;
sigmaerror = 1e-3;
fid = fopen([label '.log'], 'w');
rank = 30;

% Output file names (with the .rsf extension)
fname_jitdata = [label '_jitdata.rsf'];
fname_adjrecov = [label '_adjrecov.rsf'];
fname_L1recov = [label 'recov.rsf'];
fname_L1diff = [label 'diff.rsf'];

% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

