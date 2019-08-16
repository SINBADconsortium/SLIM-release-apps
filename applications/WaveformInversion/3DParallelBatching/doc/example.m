% Which fwi_exp##.m to display, corresponding to the same model
% fwi_exp1 - full data experiment
% fwi_exp2 - stochastic data experiment
fwi_experiments = [1 2];

base_dir = pwd; base_dir = base_dir(1:end-length('doc'));
save_dir = [base_dir 'results/'];
model_dir = [base_dir 'models/'];
fig_dir = [];

% Slices + depths to display
plots = {'z',500, 'z',1000,'x',1500,'y',1000};

% Plot figures
plot_fwi_exp(fwi_experiments,save_dir,model_dir,fig_dir,plots);