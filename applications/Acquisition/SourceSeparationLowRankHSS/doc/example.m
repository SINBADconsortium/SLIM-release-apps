%% Source separation via SVD-free rank minimization: examples and results
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: June, 2014

%%

% Set paths
curdir = pwd;
basedir = curdir(1:end-4);
datadir = [basedir '/data'];
resultsdir = [basedir '/results'];

% Load previously computed results
load([resultsdir '/SourceSep_jit1src_HSS_params.mat']) 

% Plotting parameters
dtype = 'time';
index.r = 60;
index.s = 5;
cax.freq = 5e3;
cax.time = 1e2;
cmap = 'gray';

%%

%% Original data (common shot gather)

% Source 1
makefigure([datadir '/data_zs10m.su'], dtype, index.s, cax.time, cmap, nt, dt, nr, ns)

% Source 2 (time-delayed)
makefigure([resultsdir '/' fname_D2_shift], dtype, index.s, cax.time, cmap, nt, dt, nr, ns)

% Blended shot gather 
makefigure([resultsdir '/' fname_D_blend_TRS], dtype, index.s, cax.time, cmap, nt, dt, nr, ns)

%%

%% Source separation

% Source 1
makefigure([resultsdir '/' fname_D1recov_TRS], dtype, index.s, cax.time, cmap, nt, dt, nr, ns)

% Source 2
makefigure([resultsdir '/' fname_D2recov_TRS], dtype, index.s, cax.time, cmap, nt, dt, nr, ns)

