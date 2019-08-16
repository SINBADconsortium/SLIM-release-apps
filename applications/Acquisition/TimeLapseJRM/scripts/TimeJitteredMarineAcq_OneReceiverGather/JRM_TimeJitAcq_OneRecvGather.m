% This script simulates the time-jittered (blended) marine acquisition scenario, and recovers conventional (deblended) data by sparse inversion (via one-norm minimization).
%
% Acquisition scenario: time-jittered OBC marine acquisition with one source vessel and two airgun arrays (i.e., airgun arrays fire at jittered instances in time which translate to jittered shot locations for a given speed of the source vessel)
%
% Overlap in the acquisition matrices: 20%
%
% The results are stored in the path defined in the script setpath.m


% Set the parameters (by running the JRM_TimeJitAcq_params.m script)
% NOTE: all the parameters are saved in a .mat file in the corresponding expdir
JRM_TimeJitAcq_params;

% Load data
load(fname_input);
D_base = D1;
D_mon = D2;

clear D1 D2

% Select a subset
subset = ind_start : (ind_start + ns - 1);
D1 = squeeze(D_base(:,subset,subset));
D2 = squeeze(D_mon(:,subset,subset));

% Select a receiver gather (from the subset)
D1 = squeeze(D1(:,recv_ind,:));
D2 = squeeze(D2(:,recv_ind,:));

clear D_base D_mon

% Original 4D signal
signal_4D = D1 - D2;


% Parameters for jittered acquisition
[jitacq1, jitacq2] = jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams);

% Dimensions for jittered acquisition 
nr = 1;
jitdim1 = [(length(jitacq1.tfirejitgrid) + nt - 1)  nr];
jitdim2 = [(length(jitacq2.tfirejitgrid) + nt - 1)  nr];

% Angular frequency samples
wn = [0:nt/2,-nt/2+1:-1]*2*pi/(nt*dt); 

% Sampling operator 
RM1 = opJittshift([nt nr ns], jitdim1, jitacq1, wn');
RM2 = opJittshift([nt nr ns], jitdim2, jitacq2, wn');


% Generate jittered acquisition measurements
y1 = RM1*D1(:);
y2 = RM2*D2(:);
write_ACQdata(fname_base_jitdata, reshape(y1, jitdim1), dt, dr, ds)
write_ACQdata(fname_mon_jitdata, reshape(y2, jitdim2), dt, dr, ds)


% Applying the adjoint of the sampling operator
D1est_adj = reshape((RM1')*y1, [nt ns]);
D2est_adj = reshape((RM2')*y2, [nt ns]); 
signal_4D_adj = D1est_adj - D2est_adj;

% Save adjoint-applied data
rsf_write_all(fname_base_adjrecov, {'out=stdout'}, D1est_adj, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})
rsf_write_all(fname_mon_adjrecov, {'out=stdout'}, D2est_adj, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})
rsf_write_all(fname_4Dsignal_adjrecov, {'out=stdout'}, signal_4D_adj, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})


% Operators for joint recovery 

% opCurvelet(M, N, nbscales, nbangles, finest, ttype, is_real)
C = opCurvelet(nt, ns, max(1,ceil(log2(min(nt,ns)) - 3)), 16, 1, 'ME', 0);

% Measurement operator
A1 = RM1*C'; 
A2 = RM2*C'; 
A = [opStack(A1,A2) opBlockDiag(A1,A2)];

% Recovery operator
S = [opStack(C',C') opBlockDiag(C',C')];


% Solve the one-norm recovery problem
% spgl1 : (A, b, tau, sigma, x, options)
% xest : estimated (synthesis) curvelet coefficients
xest = spgl1(A, [y1;y2], spgl1_tau, spgl1_sigma, spgl1_x, options);

% Recover data
Dest = S*xest;
D1est = reshape(Dest(1:nt*nr*ns), [nt ns]); 
D2est = reshape(Dest(nt*nr*ns+1:end), [nt ns]);
signal_4D_recov = D1est - D2est;
diff_D1 = D1 - D1est;
diff_D2 = D2 - D2est;
diff_signal_4D = signal_4D - signal_4D_recov;

% Save estimated curvelet coefficients
rsf_write_all(fname_xest, {'out=stdout'}, xest)

% Save recovered data
rsf_write_all(fname_base_L1recov, {'out=stdout'}, D1est, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})
rsf_write_all(fname_mon_L1recov, {'out=stdout'}, D2est, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})
rsf_write_all(fname_4Dsignal_L1recov, {'out=stdout'}, signal_4D_recov, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})

% Save difference data
rsf_write_all(fname_base_L1diff, {'out=stdout'}, diff_D1, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})
rsf_write_all(fname_mon_L1diff, {'out=stdout'}, diff_D2, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})
rsf_write_all(fname_4Dsignal_L1diff, {'out=stdout'}, diff_signal_4D, [nt ns], [0 0], {'Time' 'Source position'}, {'s' 'm'})


% Return to the scripts directory
cd(curdir)

