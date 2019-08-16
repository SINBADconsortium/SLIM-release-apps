%%
% Written by Curt Da Silva (curtd@math.ubc.ca), July 2015.

%% 3D FWI 
%
% This script displays the results of a very simple 3D FWI inversion
%
% The modeling used in this example is described in
% <https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/>.
%
% *System requirements:* 
%
% * This script was tested using Matlab 2014b with the parallel computing
% toolbox.
%
% * 

%% Models
cur_dir = pwd; 
base_dir = [cur_dir '/'];
results_dir = [base_dir 'results/'];

load([results_dir 'edam_fwi.mat']);

plot_slices = @(x) slice3D(x,model,1500,1500,[]);

plot_slices(v); ax = caxis;
title('true v');  colorbar;

plot_slices(v0); caxis(ax);
title('initial v');  colorbar;

plot_slices(vest); caxis(ax);
title('inverted v');  colorbar;

figure; imagesc(real(init_res)); title('Initial data residual'); ax = caxis;

figure; imagesc(real(final_res)); title('Final data residual'); caxis(ax);