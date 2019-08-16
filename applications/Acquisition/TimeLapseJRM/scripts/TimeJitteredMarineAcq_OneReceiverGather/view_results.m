% This script sets the parameters to view the results of the JRM_TimeJitAcq_OneRecvGather.m script.

% Use the experiment label
label = 'JRM_TimeJitAcq_1boat';

% Set paths
setpath;
expdir = [resultsdir '/' label];
if ~exist(expdir,'dir')
  error('The corresponding experiment directory does not exist. Please check if you have run the main script (jitacq_deblending.m)');
end
cd(expdir);

% Seismic data parameters
% (nt, nr, ns): time, receiver, and source samples 
nt = 512;
nr = 100; 
ns = 100; 

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
% see function: jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams)
% NOTE: set fig = 'no' when running the script in batch mode
rndfactor = [1000 100];
p = 4;
rseed = [1227 6782 2675 4387];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = [10.0 10.0 0.0];
delayboat = 0;

% File names (with .rsf extension, same as defined in JRM_TimeJitAcq_params.m)
fname_base_adjrecov = [label '_base_adjrecov.rsf'];
fname_mon_adjrecov = [label '_mon_adjrecov.rsf'];
fname_4Dsignal_adjrecov = [label '_4Dsignal_adjrecov.rsf'];
fname_base_L1recov = [label '_base_L1recov.rsf'];
fname_mon_L1recov = [label '_mon_L1recov.rsf'];
fname_4Dsignal_L1recov = [label '_4Dsignal_L1recov.rsf'];
fname_base_L1diff = [label '_base_L1diff.rsf'];
fname_mon_L1diff = [label '_mon_L1diff.rsf'];
fname_4Dsignal_L1diff = [label '_4Dsignal_L1diff.rsf'];

% Plotting parameters
caxdata = 0.5;
cax4Dsignal = 0.05;
cmap = 'gray';

%====================================================================================================================%

% View results
% You do NOT need to make any changes in this section

% Conventional acquisition
convacq_2arrays(flipflop, tfireint_min, ns, dt, boatspeed_conv, fig, []);

% Time-jittered acquisition
jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, []);

% x-axis label for receiver gathers
xlab = 'Source position (m)';
xpos = 0 : ds : (ns-1)*ds;
x = 1:ns;
axis_gather = 0 : 250 : 1000;
for k = 1:length(axis_gather); xtick(k) = x(xpos == axis_gather(k)); end  
xticklabel = axis_gather;

% y-axis label for receiver gathers
ylab = 'Time (s)';
ytick = 125:125:nt;
yticklabel = ytick*dt;


% Data recovered after applying the adjoint of the sampling operator
figure; imagesc(real(rsf_read_all(fname_base_adjrecov)), [-1 1]*caxdata); colormap(cmap); 
title('Baseline (adjoint recovery)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

figure; imagesc(real(rsf_read_all(fname_mon_adjrecov)), [-1 1]*caxdata); colormap(cmap); 
title('Monitor (adjoint recovery)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

figure; imagesc(real(rsf_read_all(fname_4Dsignal_adjrecov)), [-1 1]*cax4Dsignal); colormap(cmap); 
title('4-D signal (adjoint recovery)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)


% Data recovered after sparse inversion (via L1 minimization)
figure; imagesc(real(rsf_read_all(fname_base_L1recov)), [-1 1]*caxdata); colormap(cmap); 
title('Baseline (L1 recovery)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

figure; imagesc(real(rsf_read_all(fname_mon_L1recov)), [-1 1]*caxdata); colormap(cmap); 
title('Monitor (L1 recovery)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

figure; imagesc(real(rsf_read_all(fname_4Dsignal_L1recov)), [-1 1]*cax4Dsignal); colormap(cmap); 
title('4-D signal (L1 recovery)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)


% Difference between original and L1-recovered data 
figure; imagesc(real(rsf_read_all(fname_base_L1diff)), [-1 1]*caxdata); colormap(cmap); 
title('Baseline (L1 difference)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

figure; imagesc(real(rsf_read_all(fname_mon_L1diff)), [-1 1]*caxdata); colormap(cmap); 
title('Monitor (L1 difference)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

figure; imagesc(real(rsf_read_all(fname_4Dsignal_L1diff)), [-1 1]*cax4Dsignal); colormap(cmap); 
title('4-D signal (L1 difference)'); xlabel(xlab); ylabel(ylab);
set(gca, 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

