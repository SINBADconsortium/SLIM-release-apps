%% Rank minimization based source-separation in time-jittered 2D ocean-bottom marine acquisition: examples and results
%
% Author: Rajiv Kumar (rakumar@eos.ubc.ca)
%
% Date: January, 2015


%% Time-jittered acquisition with 1 boat
%  
% See the scripts under |examples|.

%%

%%
% Load the parameters file
load([resultsdir '/TimeJitAcq_1boat_2array_LR/TimeJitAcq_1boat_2array_LR_params.mat']) 

%%

%%
% time-jittered acquisition scenario

fig = 'yes';
jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig);

%%

%% index along which data should be display
frame.t = 300;
frame.r = 60;
frame.s = 92;
caxmin = 100;
caxmax = 100;
cmap   = 'gray';

%%

%%
% Time-jittered (or blended) data volume:
% view 30 seconds of the jittered data volume
t1 = 130.0;
t2 = 160.0;
imageJitteredData([resultsdir '/TimeJitAcq_1boat_2array_LR/TimeJitAcq_1boat_2array_LR_jitdata.rsf'], t1, t2, dt, nr, dr, caxmin, caxmax, cmap);

%%

%%
% Recovery by conventional processing:
% apply the adjoint of the sampling operator
image2DTRSframe([resultsdir '/TimeJitAcq_1boat_2array_LR/TimeJitAcq_1boat_2array_LR_adjrecov.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%
% NOTE: an empty shot gather image implies that none of the airguns fired at that location.
%%

%% 
% Recovery by rank-minimization
image2DTRSframe([resultsdir '/TimeJitAcq_1boat_2array_LR/TimeJitAcq_1boat_2array_LRrecov.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%
% Difference
image2DTRSframe([resultsdir '/TimeJitAcq_1boat_2array_LR/TimeJitAcq_1boat_2array_LRdiff.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%% Running the code on your own data
%
% Instructions are included in the README.md file under
% |/applications/Acquisition/2DTimeJitteredOBS-LR/example|.

