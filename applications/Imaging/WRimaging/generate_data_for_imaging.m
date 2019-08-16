%% Generate data for Wavefield Reconstruction Imaging examples
%this is a synthetic example of WRI using the true source wavelet. Sources
%are close to the water surface and receivers are at the water bottom. Full
%source and reveir offset is used.
%
%
% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: July, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)
%
%This script will generate data for the Wavefield Reconstruction Imaging
%example 'WRImaging_example.m'. It runs for about 10 minutes on 5 nodes
%with 5 processes per node.

close all

startup
addpath([slimapps, '/applications/WaveformInversion/2DWRI/mbin']);

%% Main preferences (all user input is defined in this section)

%logistics
save_data           = 1;
save_workspace      = 0;
save_model_iters    = 1; %save model estimate at every iteration
save_model_batch    = 1; %save the model estimate at the end of every frequency batch

%data directory
data_dir = '/data/';

% Gaussian noise
add_noise = 0;
NSR   = 0.075; %noise to signal ration (2norm) ; NSR = norm(noise)/norm(signal)

%use simultaneous sources (yes=1). Redrawing means new random source
%weights are selected at the end of every iteration. Redrawing is mostly suitable for
%gradient-descent or Gauss-Newton type algorithms, but not for quasi-Newton
%algorithms.
params.sim_src        = 0;
params.n_sim_src      =20;
params.redraw_sources = 0;

params.C        = 1; %this is a source/receiver mask in case of marine-aquisition, otherwise it is 1 (for land aquisition or ocean-bottom node)
params.srcest   = 0; %use source function estimation

%define frequency batches to be used (one at a time)
nbatch=24;        
for i=1:nbatch
    f_a{i} = linspace(5+(i-1),6+(i-1),5);
end
model.freq  = unique(cell2mat(f_a));
model.nfreq = length(model.freq); nfreq=length(model.freq); 

% Ricker wavelet (use f0=0 and t0=0 for an impulse)
model.f0 = 0;      %preak frequency
model.t0 = 0;    %time shift

%discretization options
params.nppw    = 6;   %number of grid points per wavelength (should be [6 - 12])
params.nwl_pml = 1; %physical size of the pml in terms of the longest wavelength in the start model

%computational options
params.nthreads = 4; %threads per spmd block (this is the lowest level of paralellism, used for the factorization)

% source/receiver grid and locations
%make sure these coordinates (in meters) are actually inside the
%computational domain!
model.zsrc = 25;
model.xsrc = 25:25:5725;
model.zrec = 200 ;
model.xrec = 25:25:5725;

%% startup
startup

current_dir = pwd;
addpath([current_dir,'/mbin']);
cd([current_dir,data_dir])

%% load and resample true velocity model in (m/s)
v = load(strcat(current_dir,'/data/cv.mat')); v=v.Data;

%use only a piece of the model
%v = v(1:50,1550:1600);
v = v(1:end,1535:1765);
v_true = v;
%----------
% Create background velocity model used as start model


model.d = [6 25]; %z and x axis increment in the original model
model.o = [0 0];
model.n = [size(v,1) size(v,2)];
[model.z, model.x]=odn2grid(model.o,model.d,model.n);   %original grid
z_orig = model.z;
x_orig = model.x;
units = 'm/s';
[model,int, ~,Q_scale_factor]=create_grid_params(v,model,model.freq,params,units);
v = int*v(:);
v=reshape(v,model.n);

% gridded slowness-squared values [km^2/s^2]
m = 1e6./v(:).^2;

%% Set up source and receiver arrays

model.nrec = length(model.xrec);
model.nsrc = length(model.xsrc);

% source matrix
Q = speye(length(model.xsrc));

%% Set up the generation of synthetic data in the true velocity model

%the lines below can be used to generate data and to add noise
 D           = F(m,Q_scale_factor.*Q,model,params);
 if save_data==1; Ds=gather(D); save('D_6ppw_5_29Hz_fulloffset_imaging.mat','Ds'); clear Ds;  end;
 
%add noise?
if add_noise==1;
    
    noise = randn(size(D,1),1)+1i*randn(size(D,1),1);
    noise = noise./norm(noise);
    noise = noise*NSR*norm(D(:));
    D     = D+noise;

if save_data==1; Ds=gather(D); save('D_noisy','Ds'); clear Ds;    end;
end
