% Parameters template
% This script is a template to set the parameters for time-jittered (blended) marine acquisition, and 
% conventional data recovery by sparse inversion (via one-norm minimization). 
% The script contains dummy names/values for all the variables. 
% Replace the dummy names/values with your own. Please contact us if you have any questions.


% Set a label for the experiment
label = 'test';


% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);


% Input data file
% NOTE: data file can be in SEGY/SU/RSF/MAT format;
%       simultaneously, make changes in the jitacq_deblending_template.m 
%       script to use the corresponding function in MATLAB to read the data
fname_input = 'testdata.rsf';  


% Seismic data parameters
% (nt_org, nr_org, ns_org): time, receiver, and source samples for all the data
%                           [not required if all the data is being used]
% (nt, nr, ns): time, receiver, and source samples for a subset of the data
% NOTE: nt must be a power of 2 (required by the FFT operation performed within the opSplineWaveletSPOT operator) 
nt_org = 1024;
nr_org = 178;
ns_org = 178;
nt = 512;
nr = 128; 
ns = 128; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 25.0;     
ds = 25.0;     


% Parameters for time-jittered acquisition
% see function: jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig)
% NOTE: set fig = 'no' when running the script in batch mode
p = 2;
nboats = 1;
randseed = [3 1];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = 10.0;
delayboat = 0;
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
fname_jitdata = 'jitdata.rsf';
fname_adjrecov = 'adjrecov.rsf';
fname_L1recov = 'L1recov.rsf';
fname_L1diff = 'L1diff.rsf';


% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

