%% Joint recovery method: examples and results 
%
% Author: Haneet Wason (hwason@eos.ubc.ca)
%
% Date: June, 2014

%%

% Set paths
curdir = pwd;
basedir = curdir(1:end-4);
datadir = [basedir '/data/TimeJitteredMarineAcq'];
resultsdir = [basedir '/results/TimeJitteredMarineAcq_OneReceiverGather'];

% Load previously computed results
load([resultsdir '/JRM_TimeJitAcq_1boat_params.mat']) 

% Plotting parameters
caxdata = 0.5;
cax4Dsignal = 0.05;
cmap = 'gray';

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

%%

%% Conventional vs. jittered (blended) acquisition

% Conventional acquisition
flipflop = 'yes';
tfireint_min = 10.0;
boatspeed_conv = ds/tfireint_min;
fig = 'yes';
convacq_2arrays(flipflop, tfireint_min, ns, dt, boatspeed_conv, fig, []);

% Time-jittered acquisition for baseline and monitor surveys
jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, []);


%%

%% Original data (common receiver gather)

% Load data
load('../data/TimeJitteredMarineAcq/data_4D.mat');
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

% Baseline
figure; imagesc(D1, [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('True baseline');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% Monitor
figure; imagesc(D2, [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('True monitor');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% 4-D signal
figure; imagesc(signal_4D, [-1 1]*cax4Dsignal); colormap(cmap); xlabel(xlab); ylabel(ylab); title('True 4-D signal');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

%%

%% Apply the adjoint of the sampling operator to blended data

% Baseline
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_base_adjrecov])), [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('Baseline');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% Monitor
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_mon_adjrecov])), [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('Monitor');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% 4-D signal
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_4Dsignal_adjrecov])), [-1 1]*cax4Dsignal); colormap(cmap); xlabel(xlab); ylabel(ylab); title('4-D signal');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

%%

%% Recovered data

% Baseline
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_base_L1recov])), [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('Baseline');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% Monitor
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_mon_L1recov])), [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('Monitor');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% 4-D signal
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_4Dsignal_L1recov])), [-1 1]*cax4Dsignal); colormap(cmap); xlabel(xlab); ylabel(ylab); title('4-D signal');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

%%

%% Residual : Original - Recovered

% Baseline
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_base_L1diff])), [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('Baseline');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% Monitor
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_mon_L1diff])), [-1 1]*caxdata); colormap(cmap); xlabel(xlab); ylabel(ylab); title('Monitor');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

% 4-D signal
figure; imagesc(real(rsf_read_all([resultsdir '/' fname_4Dsignal_L1diff])), [-1 1]*cax4Dsignal); colormap(cmap); xlabel(xlab); ylabel(ylab); title('4-D signal');
set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Xtick', xtick, 'XTickLabel', xticklabel, 'Ytick', ytick, 'YTickLabel', yticklabel)

