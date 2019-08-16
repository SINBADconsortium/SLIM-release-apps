% This script sets the parameters for time-jittered (blended) marine acquisition, and conventional data
% recovery by sparse inversion (via one-norm minimization). 
%
% Acquisition scenario: time-jittered OBC marine acquisition with two source vessels and two airgun arrays each


% Set a label for the experiment
label = 'TimeJitAcq_2boats';


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
% (nt, nr, ns): time, receiver, and source samples for a subset of the data
% NOTE: for 'SuezShots125.355shots.su', nt_org = 1024, nr_org = ns_org = 355, dt = 0.004, dr = ds = 12.5
%       nt must be a power of 2 (required by the FFT operation performed within the opSplineWaveletSPOT operator) 
nt_org = 1024;
nr_org = 355;
ns_org = 355;
nt = 1024;
nr = 128; 
ns = 128; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004; 
dr = 12.5;     
ds = 12.5;      


% Parameters for time-jittered acquisition
% see function: jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig)
% NOTE: set fig = 'no' when running the script in batch mode
p = 4;
nboats = 2;
randseed = [8 5 3 7];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = 10.0;
delayboat = 300.0;
fig = 'no';


% Flag to perform dottest for the sampling operator
test = 0;


% Set options for the SPGL1 solver
% see solver: spgl1(A, b, tau, sigma, x, options)
spgl1_tau = 0;
spgl1_sigma = 0;
spgl1_x = [];
options.fid = fopen([label '.log'], 'w');
options.verbosity = 1;
options.iterations = 200;
options.optTol = 1e-04;


% Output file names (with the .rsf extension)
fname_jitdata = [label '_jitdata.rsf'];
fname_adjrecov = [label '_adjrecov.rsf'];
fname_L1recov = [label '_L1recov.rsf'];
fname_L1diff = [label '_L1diff.rsf'];


% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

