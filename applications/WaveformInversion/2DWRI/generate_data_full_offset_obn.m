%% Generate data for the examples in the WRI application. - BG compass model
% the same data can also be obtained by running the SConstruct script in
% the data folder.
% Sources are close to the water surface and receivers are at the water bottom. Full
% source and reveir offset is used. This example will run about <5 minutes
% on 1 node with 5 processes per node on ~2015 hardware).

% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: July, 2015, last updated: February 2016

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)

%This script will generate data and invert using the same modeling code.
%Model is the Compass model (part of it). Correct source function
%is used.

close all

%% Main preferences (all user input is defined in this section)

%logistics
save_data           = 1;

%data folder, files and results directory
data_dir = '/data/';

%define frequency batches to be used
nbatch=24;
for i=1:nbatch
    f_a{i} = linspace(5+(i-1),6+(i-1),4);
end
model.freq  = unique(cell2mat(f_a));
model.nfreq = length(model.freq); nfreq=length(model.freq);

% Ricker wavelet (use f0=0 and t0=0 for an impulse)
model.f0 = 15;      %preak frequency
model.t0 = 0;       %time shift

% source/receiver grid and locations
%make sure these coordinates (in meters) are actually inside the
%computational domain!
model.zsrc = 25;
model.xsrc = 25:50:5725;
model.zrec = 200 ;
model.xrec = 25:50:5725;

%computational options
%params.nthreads = 1; %threads per spmd block (this is the lowest level of paralellism, used for the LU factorization)
params.C        = 1; %this is a source/receiver mask in case of marine-aquisition, otherwise it is 1 (for land aquisition or ocean-bottom node)

if isfield(params,'pdefunopts')==0
    params.pdefunopts = PDEopts2D();
end
if isfield(params,'lsopts')==0
    solve_opts = LinSolveOpts();
    solve_opts.solver = LinSolveOpts.SOLVE_LU;
    params.lsopts = solve_opts;
end
    
%% startup
startup

current_dir = pwd;
addpath([current_dir,'/mbin']);
cd([current_dir,data_dir])

%% load and resample true velocity model in (m/s)
v = load(strcat(current_dir,'/data/cv.mat')); v=v.Data;

%use only a piece of the model
v = v(1:end,1535:1765);

model.d = [6 25]; %z and x axis increment (in m) in the original model (6 & 25 is for the BG compass)
model.o = [0 0];
model.n = [size(v,1) size(v,2)];
[model.z, model.x]=odn2grid(model.o,model.d,model.n);   %original grid

% v=reshape(v,model.n);
% plot velocity model to be used
%figure;
%imagesc(model.x,model.z,v);axis equal tight;colorbar;
%xlabel('x [m]');ylabel('z [m]');title('True velocity model [m/s]');

% gridded slowness-squared values [km^2/s^2]
m = 1e6./v(:).^2;

%% Set up the generation of synthetic data in the true velocity model
model.n(3)=1;
model.d(3)=1;
model.o(3)=0;
model.nrec = length(model.xrec);
model.nsrc = length(model.xsrc);

% source matrix
Q = speye(length(model.xsrc));

% generate data
D           = F(m,Q,model,params);
Ds=gather(D);


%save
datafilename=['BGcompass_obn_fulloffset_',num2str(min(model.freq)),'_',num2str(max(model.freq)),'Hz'];
data_info.freq=model.freq;
data_info.xrec=model.xrec;
data_info.zrec=model.zrec;
data_info.xsrc=model.xsrc;
data_info.zsrc=model.zsrc;
data_info.nrec=model.nrec;
data_info.nsrc=model.nsrc;
data_info.t0=model.t0;
data_info.f0=model.f0;
save(datafilename,'Ds','data_info');
