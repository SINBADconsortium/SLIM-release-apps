% This script sets the parameters for over/under blended marine acquisition and deblending (or source separation) by promoting low rank via nuclear-norm minimization.

% Set a label for the experiment
label = 'SourceSep_jit1src_HSS';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% Input data file
fname_input1 = [datadir1 '/data_zs10m.su'];
fname_input2 = [datadir2 '/data_zs15m.su'];

% Seismic data parameters
% (nt, nr, ns): number of time, receiver, and source samples
nt = 1024;
nr = 129;
ns = 129;

% (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;
dfreq = 1/(nt*dt);
dr = 25.0;
ds = 25.0;

% Use only one side of the FX spectra, since the other half is symmetric
if mod(nt,2) == 0 
   nfreq = floor((nt/2)+1);
else
   nfreq = ceil(nt/2);
end

% Set the random seed
% see function: gen_blended_data
rseed = 1;

% Set the number of HSS levels
HSS_level = 1;

% Set options for the SPGL1 solver
opts = spgSetParms('verbosity', 1, ...
                   'decTol', 1e-6, ...
                   'project', @TraceNorm_project_macq, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual_macq, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 1000, ...
                   'weights', []);
opts.funPenalty = @funLS;
params.mode = 2;
params.depth = 5;
params.watervel = 1480;
params.ls = 1 ;

% Frequency axis
freq_axis = dfreq*(0:nfreq-1);

% Range of frequencies with significant energy
freq_min = 0;
freq_max = 80;

% Set minimum and maximum rank values 
% (for the lowest and highest frequency slices, respectively)
rk_min = 10;
rk_max = 50;

% Output file names (with the .rsf extension)
fname_D1 = [label '_D1.rsf'];
fname_D2_shift = [label '_D2_shifted.rsf'];
fname_D_blend_FRS = [label '_blended_FRS.rsf'];
fname_D_blend_TRS = [label '_blended_TRS.rsf'];
fname_D1recov_FRS = [label '_src1recov_FRS.rsf'];
fname_D1recov_TRS = [label '_src1recov_TRS.rsf'];
fname_D2recov_FRS = [label '_src2recov_FRS.rsf'];
fname_D2recov_TRS = [label '_src2recov_TRS.rsf'];

% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

