% This script sets the parameters to view the results of the run_TimeJitAcq_Deblending.m script.

% Use the experiment label
label = 'TimeJitAcq_OffTheGrid';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
   error('The corresponding experiment directory does not exist. Please check if you have run the main script (run_TimeJitAcq_Deblending.m)');
end
cd(expdir);

% Seismic data parameters
% (nt, nr, ns): time, receiver, and source samples 
nt = 1024;
nr = 128; 
ns = 128; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 12.5;     
ds = 12.5;     

% Parameters for conventional acquisition
flipflop = 'yes';
tfireint_min = 10.0;
boatspeed_conv = ds/tfireint_min;
fig = 'yes';

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
fig = 'yes';
figparams = [];

% File names (with .rsf extension)
fname_irregdata = [label '_irregdata.rsf'];
fname_jitdata = [label '_jitdata.rsf'];
fname_adjrecov = [label '_adjrecov.rsf'];
fname_L1recov = [label '_L1recov.rsf'];
fname_L1diff = [label '_L1diff.rsf'];

% Parameters to view a section of the jittered data volume
% see function: imageJitteredData(fname, t1, t2, dt, nr, dr, caxmin, caxmax, cmap)
t1 = 100.0;
t2 = 130.0;
caxmin = 50;
caxmax = 50;
cmap = 'gray';

%====================================================================================================================%

% View results
% You do NOT need to make any changes in this section

% Conventional acquisition
convacq_2arrays(flipflop, tfireint_min, ns, dt, boatspeed_conv, fig, []);

% Time-jittered acquisition (off the grid)
jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, [])

% A time, receiver, and source gather of the irregular data
image2DTRSframe(fname_irregdata, frame, nt, dt, caxmin, caxmax, cmap);

% A section of the jittered data volume
imageJitteredData(fname_jitdata, t1, t2, dt, nr, dr, caxmin, caxmax, cmap); 

% A time, receiver, and source gather of the data recovered after applying the adjoint of the sampling operator
image2DTRSframe(fname_adjrecov, frame, nt, dt, caxmin, caxmax, cmap);

% A time, receiver, and source gather of the data recovered after sparse inversion (via L1 minimization)
image2DTRSframe(fname_L1recov, frame, nt, dt, caxmin, caxmax, cmap);

% A time, receiver, and source gather of the difference between original and L1-recovered data 
image2DTRSframe(fname_L1diff, frame, nt, dt, caxmin, caxmax, cmap);

