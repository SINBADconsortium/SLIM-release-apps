% This script sets the parameters for time-jittered (blended) marine acquisition, and conventional data recovery (or deblending) by sparse inversion (via one-norm minimization) on non-equispaced spatial grids 
%
% Acquisition scenario --- time-jittered OBC marine acquisition ('off the grid') with one source vessel and two airgun arrays 


% Set a label for the experiment
label = 'TimeJitAcq_1boat_50to12pt5m_NFDCT';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% Input data file
fname_input = [datadir '/SuezShots125.355shots.su'];  

% Seismic data parameters
% (nt_org, nr_org, ns_org): time, receiver, and source samples for all the data
nt_org = 1024;
nr_org = 355;
ns_org = 355;

% Select a subset of the data to work on
% (nt, nr, ns): time, receiver, and source samples for a subset of the data
% NOTE - nt must be a power of 2 (required by the FFT operation performed within the opSplineWaveletSPOT operator) 
nt = 1024;
nr = 128; 
ns = 128; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 12.5;     
ds = 12.5;     

% Parameters for time-jittered acquisition
% see function: jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams)
% NOTE: set fig = 'no' when running the script in batch mode
rndfactor = [1000 100];
p = 4;
nboats = 1;
rseed = [68 71]; 
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = 10.0;
delayboat = 0;
fig = 'no';
figparams = [];

% Flag to perform dottest for the sampling operator
% NOTE - 0 ==> no, 1 ==> yes
test = 0;

% Set options for the SPGL1 solver
% see function: spgl1(A, b, tau, sigma, x, options)
opt.spgl1_tau = 0;
opt.spgl1_sigma = 0;
opt.spgl1_x = [];
options.fid = fopen([label '.log'], 'w');
options.verbosity = 1;
options.iterations = 200;
options.optTol = 1e-04;
options.ignorePErr = 1;

% Output file names
fname_irregdata = [label '_irregdata.rsf'];
fname_jitdata = [label '_jitdata.rsf'];
fname_adjrecov = [label '_adjrecov.rsf'];
fname_L1recov = [label '_L1recov.rsf'];
fname_L1diff = [label '_L1diff.rsf'];

% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

