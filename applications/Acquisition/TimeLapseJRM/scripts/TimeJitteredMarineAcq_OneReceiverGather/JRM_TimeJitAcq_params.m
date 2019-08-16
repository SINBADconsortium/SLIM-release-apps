% This script sets the parameters for time-jittered (blended) marine acquisition, and conventional data recovery (or deblending) by sparse inversion (via one-norm minimization) 
%
% Acquisition scenario: time-jittered OBC marine acquisition with one source vessel and two airgun arrays


% Set a label for the experiment
label = 'JRM_TimeJitAcq_1boat';


% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);


% Input data file
fname_input = [datadir '/data_4D.mat'];  


% Seismic data parameters
% (nt_org, nr_org, ns_org): time, receiver, and source samples for all the data
% nt must be a power of 2 (required by the FFT operation performed within the opSplineWaveletSPOT operator) 
nt_org = 512;
nr_org = 151;
ns_org = 151;

% Select a subset of the data
% recv_ind: receiver index (to work on a receiver gather)
% ind_start: starting shot index (work with 'ns' shots defined below)
recv_ind = 50;
ind_start = 30;

% (nt, nr, ns): time, receiver, and source samples for a subset of the data
nt = 512;
nr = 100; 
ns = 100; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 12.5;     
ds = 12.5;     


% Parameters for time-jittered acquisition
% see function: jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams)
% NOTE: set fig = 'no' when running the script in batch mode
rndfactor = [1000 100];
p = 4;
rseed = [1227 6782 2675 4387];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = [10.0 10.0 0.0];
delayboat = 0;
fig = 'no';
figparams = [];


% Flag to perform dottest for the sampling operator
test = 0;


% Set options for the SPGL1 solver
% see function: spgl1(A, b, tau, sigma, x, options)
spgl1_tau = 0;
spgl1_sigma = 0;
spgl1_x = [];
options.fid = fopen([label '.log'], 'w');
options.verbosity = 1;
options.iterations = 500;
options.optTol = 1e-04;
options.ignorePErr = 1;


% Output file names (with the .rsf extension)
fname_base_jitdata = [label '_base_jitdata.rsf'];
fname_mon_jitdata = [label '_mon_jitdata.rsf'];
fname_base_adjrecov = [label '_base_adjrecov.rsf'];
fname_mon_adjrecov = [label '_mon_adjrecov.rsf'];
fname_4Dsignal_adjrecov = [label '_4Dsignal_adjrecov.rsf'];
fname_xest = [label '_xest.rsf'];
fname_base_L1recov = [label '_base_L1recov.rsf'];
fname_mon_L1recov = [label '_mon_L1recov.rsf'];
fname_4Dsignal_L1recov = [label '_4Dsignal_L1recov.rsf'];
fname_base_L1diff = [label '_base_L1diff.rsf'];
fname_mon_L1diff = [label '_mon_L1diff.rsf'];
fname_4Dsignal_L1diff = [label '_4Dsignal_L1diff.rsf'];


% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

