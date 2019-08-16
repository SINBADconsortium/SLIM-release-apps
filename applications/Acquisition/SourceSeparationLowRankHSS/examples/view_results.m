% This script sets the parameters to view the results of the SourceSep_HSS.m script.

% Use the experiment label
label = 'SourceSep_jit1src_HSS';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
  error('The corresponding experiment directory does not exist. Please check if you have run the main script (jitacq_deblending.m)');
end
cd(expdir);

% Seismic data parameters
% (nt, nr, ns): time, receiver, and source samples 
nt = 1024;
nr = 129; 
ns = 129; 

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dfreq = 1/(nt*dt);
dr = 25.0;     
ds = 25.0;     

% File names (with .rsf extension)
fname_D1 = [label '_D1.rsf'];
fname_D2_shift = [label '_D2_shifted.rsf'];
fname_D_blend_FRS = [label '_blended_FRS.rsf'];
fname_D_blend_TRS = [label '_blended_TRS.rsf'];
fname_D1recov_FRS = [label '_src1recov_FRS.rsf'];
fname_D1recov_TRS = [label '_src1recov_TRS.rsf'];
fname_D2recov_FRS = [label '_src2recov_FRS.rsf'];
fname_D2recov_TRS = [label '_src2recov_TRS.rsf'];

% Plotting parameters
index.f = 40;
index.t = 100;
index.r = 60;
index.s = 5;
cax.freq = 5e3;
cax.time = 1e2;
cmap = 'gray';

%====================================================================================================================%

% Plots

% Frequency domain data
dtype = 'frequency';
makefigure(fname_D_blend_FRS, dtype, index.f, cax.freq, cmap, nt, dt) % Blended data
makefigure(fname_D1recov_FRS, dtype, index.f, cax.freq, cmap, nt, dt) % Recovered source 1
makefigure(fname_D2recov_FRS, dtype, index.f, cax.freq, cmap, nt, dt) % Recovered source 2

% Time domain data
dtype = 'time';
makefigure(fname_D1, dtype, index.s, cax.time, cmap, nt, dt); % Original shot gather (source 1)
makefigure(fname_D2_shift, dtype, index.s, cax.time, cmap, nt, dt); % Original shot gather (source 2)
makefigure(fname_D_blend_TRS, dtype, index.s, cax.time, cmap, nt, dt); % Blended shot gather
makefigure(fname_D1recov_TRS, dtype, index.s, cax.time, cmap, nt, dt); % Recovered source 1
makefigure(fname_D2recov_TRS, dtype, index.s, cax.time, cmap, nt, dt); % Recovered source 2

