% This script sets the parameters for over/under blended marine acquisition and deblending (or source separation) by sparsity promotion.

% Set a label for the experiment
label = 'SourceSepL1';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% Data file names - blended data, time delays for each source, and time-shifted data
% FRS: data in frequency-receiver-source coordinates
fname_D_blend = [datadir '/blended_data.rsf'];
fname_D_blend_FRS = [datadir '/blended_data_FRS.rsf'];
fname_delay_S1 = [datadir '/tdelay_src1.rsf'];
fname_delay_S2 = [datadir '/tdelay_src2.rsf'];
fname_D1_shift = [datadir '/shifted_src1.rsf'];
fname_D2_shift = [datadir '/shifted_src2.rsf'];
fname_D1_shift_FRS = [datadir '/shifted_src1_FRS.rsf'];
fname_D2_shift_FRS = [datadir '/shifted_src2_FRS.rsf'];

% Data dimensions
% (nt, nr, ns): number of time, receiver, and source samples
nt = 1250; 
nr = 231;
ns = 231;

% Sampling intervals
% (dt, dr, ds): time, receiver, and source sampling interval
dt = 0.004;
dfreq = 1/(nt*dt);
dr = 20;
ds = 20;

% Use only one side of the FX spectra, since the other half is symmetric
if mod(nt,2) == 0
   nfreq = floor((nt/2)+1);
else
   nfreq = ceil(nt/2);
end

% Frequency axis
freq_axis = dfreq*(0:nfreq-1);

% Range of frequencies to separate blended data
% NOTE - all frequencies used for source separation 
freq_min = 0;
freq_max = freq_axis(end);

% Set options for the SPGL1 solver
% see function: spgl1(A, b, tau, sigma, x, options)
opt.spgl1_tau = 0;
opt.spgl1_sigma = 0;
opt.spgl1_x = [];
options.fid = fopen([label '.log'], 'w');
options.verbosity = 1;
options.iterations = 350;
options.optTol = 1e-04;
options.ignorePErr = 1;

% Output file names (with the .rsf extension)
% FRS: data in frequency-receiver-source coordinates
% TRS: data in time-receiver-source coordinates
fname_D1recov_FRS = [label '_src1recov_FRS.rsf'];
fname_D2recov_FRS = [label '_src2recov_FRS.rsf'];
fname_D1shiftrecov_FRS = [label '_shifted_src1recov_FRS.rsf'];
fname_D2shiftrecov_FRS = [label '_shifted_src2recov_FRS.rsf'];
%fname_D1recov_TRS = [label '_src1recov_TRS.rsf'];
%fname_D2recov_TRS = [label '_src2recov_TRS.rsf'];
%fname_D1shiftrecov_TRS = [label '_shifted_src1recov_FRS.rsf'];
%fname_D2shiftrecov_TRS = [label '_shifted_src2recov_FRS.rsf'];


% Save the parameters in a .mat file
fname_mat = [label '_params.mat'];
save(fname_mat);

