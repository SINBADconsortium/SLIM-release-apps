%% 2D ocean-bottom marine acquisition via jittered sampling: examples and results
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: April, 2013


%% Time-jittered acquisition with 1 boat
%  
% See the scripts under |examples/TimeJitAcq_1boat/|.

%%

%%
% Load the parameters file
load([resultsdir '/TimeJitAcq_1boat/TimeJitAcq_1boat_params.mat']) 

%%

%%
% Conventional vs. time-jittered acquisition scenario
flipflop = 'no';
fig = 'yes';
convacq_2arrays(flipflop, tfireint_min, ns, dt, boatspeed, fig);
jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig);

%%

%%
% Original data
frame.t = 300;
frame.r = 60;
frame.s = 92;
caxmin = 100;
caxmax = 100;
cmap   = 'gray';
image2DTRSframe([resultsdir '/TimeJitAcq_1boat/Suez_ds25m.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%

%%
% Time-jittered (or blended) data volume:
% view 30 seconds of the jittered data volume
t1 = 130.0;
t2 = 160.0;
imageJitteredData([resultsdir '/TimeJitAcq_1boat/TimeJitAcq_1boat_jitdata.rsf'], t1, t2, dt, nr, dr, caxmin, caxmax, cmap);

%%

%%
% Recovery by conventional processing:
% apply the adjoint of the sampling operator
image2DTRSframe([resultsdir '/TimeJitAcq_1boat/TimeJitAcq_1boat_adjrecov.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%
% NOTE: an empty shot gather image implies that none of the airguns fired at that location.
%
% This step is usually followed by some kind of median filtering on common receiver or CMP gathers.
% We address the challenge of deblending by a (curvelet-based) sparsity-promoting recovery technique.

%%

%% 
% Recovery by sparse inversion (via L1 minimization)
image2DTRSframe([resultsdir '/TimeJitAcq_1boat/TimeJitAcq_1boat_L1recov.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%
% Difference
image2DTRSframe([resultsdir '/TimeJitAcq_1boat/TimeJitAcq_1boat_L1diff.rsf'], frame, nt, dt, caxmin, caxmax, cmap);



%% Time-jittered acquisition with 2 boats 
%
% See the scripts under |examples/TimeJitAcq_2boats/|.

%%

%%
% Load the parameters file
load([resultsdir '/TimeJitAcq_2boats/TimeJitAcq_2boats_params.mat']) 

%%

%%
% Time-jittered acquisition scenario
fig = 'yes';
jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig);

%%

%%
% Original data
frame.t = 300;
frame.r = 70;
frame.s = 30;
image2DTRSframe([resultsdir '/TimeJitAcq_2boats/Suez_ds12pt5m.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%

%%
% Time-jittered (or blended) data volume:
% view 30 seconds of the jittered data volume
t1 = 330.0;
t2 = 360.0;
imageJitteredData([resultsdir '/TimeJitAcq_2boats/TimeJitAcq_2boats_jitdata.rsf'], t1, t2, dt, nr, dr, caxmin, caxmax, cmap);

%%

%%
% Recovery by conventional processing:
% apply the adjoint of the sampling operator
image2DTRSframe([resultsdir '/TimeJitAcq_2boats/TimeJitAcq_2boats_adjrecov.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%
% This step is usually followed by some kind of median filtering on common receiver or CMP gathers.
% We address the challenge of deblending by a (curvelet-based) sparsity-promoting recovery technique.

%%

%% 
% Recovery by sparse inversion (via L1 minimization)
image2DTRSframe([resultsdir '/TimeJitAcq_2boats/TimeJitAcq_2boats_L1recov.rsf'], frame, nt, dt, caxmin, caxmax, cmap);

%%
% Difference
image2DTRSframe([resultsdir '/TimeJitAcq_2boats/TimeJitAcq_2boats_L1diff.rsf'], frame, nt, dt, caxmin, caxmax, cmap);



%% Running the code on your own data
%
% Template scripts are included that can be adapted to run the code on
% your own data. See the |scripts| directory.

