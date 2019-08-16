% This script is a template to set the parameters and view the results of the jitacq_deblending_template.m script.
% The script contains dummy names/values for all the variables. 
% Replace the dummy names/values with your own. Please contact us if you have any questions.


% Use the experiment label
label = 'test';


% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
   error('The corresponding experiment directory does not exist. Please check if you have run the main script (jitacq_deblending_template.m)');
end
cd(expdir);


% Seismic data parameters
% (nt, nr, ns): time, receiver, and source samples 
nt = 1024;
nr = 128; 
ns = 128; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 25.0;     
ds = 25.0;     


% Parameters for time-jittered acquisition 
% see function: jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig)
% NOTE: use the same values as used in jitacq_deblending_params_template.m, and set fig = 'yes'
p = 2;
nboats = 1;
randseed = [3 1];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = 10.0;
delayboat = 0;
fig = 'yes';


% Flag to display results
% JitData: time-jittered data volume
% AdjRecovData: data recovered after applying the adjoint of the sampling operator
% L1RecovData: data recovered by sparse inversion via one-norm minimization
% L1DiffData: difference between original and L1-recovered data 
JitData = 1;       
AdjRecovData = 1;   
L1RecovData = 1;   
L1DiffData = 1;    


% File names (with .rsf extension)
fname_jitdata = 'test_jitdata.rsf';
fname_adjrecov = 'test_adjrecov.rsf';
fname_L1recov = 'test_L1recov.rsf';
fname_L1diff = 'test_L1diff.rsf';


% Parameters to view a section of the jittered data volume
% see function: imageJitteredData(fname, t1, t2, dt, nr, dr, caxmin, caxmax, cmap)
t1 = 100.0;
t2 = 140.0;
caxmin = 100;
caxmax = 100;
cmap = 'gray';


% Parametes to view a time, receiver, and source gather of a data volume
% see function: image2DTRSframe(fname, frame, nt, dt, caxmin, caxmax, cmap)
% NOTE: use this to view AdjRecovData, L1RecovData, and L1DiffData 
frame.t = 200;
frame.r = 60;
frame.s = 101;


%====================================================================================================================%


% View results

% View the acquisition scenario
if strcmp(fig, 'yes')
   jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig);
end


% View a section of the jittered data volume
if JitData, 
   imageJitteredData(fname_jitdata, t1, t2, dt, nr, dr, caxmin, caxmax, cmap); 
end


% View a time, receiver, and source gather of the data recovered after applying the adjoint of the sampling operator
if AdjRecovData,
   image2DTRSframe(fname_adjrecov, frame, nt, dt, caxmin, caxmax, cmap);
end


% View a time, receiver, and source gather of the data recovered after sparse inversion (via L1 minimization)
if L1RecovData,
   image2DTRSframe(fname_L1recov, frame, nt, dt, caxmin, caxmax, cmap);
end


% View a time, receiver, and source gather of the difference between original and L1-recovered data 
if L1DiffData,
   image2DTRSframe(fname_L1diff, frame, nt, dt, caxmin, caxmax, cmap);
end

