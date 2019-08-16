%% Wavefield Reconstruction Imaging
%
% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: April, 2014. Some small updates in July 2015.
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)
%
% This script generates an image of the medium, given an estimate of the
% velocity model. The goal is the same as in any migration algorithm, but
% the method is new. This script is based on the WRI objective function and
% is described in [1,2] and an example is shown in [3]. This is compared to
% basic Reverse-time-migration. Example uses the true source wavelet (5-29
% Hz) using sources near the water surface and receivers on the ocean
% bottom. Source and receiver spacing is 25 m.
%
% Some computational notes:
% Runtime is about 15 minutes on 5 nodes with 5 processes (factorizations)
% per node. More than 1 core per process can be used by setting
% 'params.nthreads = 1;' to a number larger than one. This example has 97
% frequencies which are all solved on the same grid.


close all

startup
addpath([slimapps, '/applications/WaveformInversion/2DWRI/mbin']);

%% Main preferences (all user input is defined in this section)

%logistics
save_figures        = 1;
save_data           = 0;

%data folder, files and results directory
resultsdir = '/results';
%datadir    = '/data/';
%datafile   = 'D_6ppw_5_29Hz_fulloffset_imaging.mat';

% Gaussian noise
add_noise    = 0;
NSR   = 0.075; %noise to signal ration (2norm) ; NSR = norm(noise)/norm(signal)

params.C        = 1; %this is a source/receiver mask in case of marine-aquisition, otherwise it is 1 (for land aquisition or ocean-bottom node)
params.srcest   = 0; %use source function estimation

%define frequencies
model.freq  = 5:0.2:24.8;%unique(cell2mat(f_a));
model.nfreq = length(model.freq);

spmd;nl=numlabs;end;
if mod(model.nfreq,nl{1})~=0
    error('number of frequencies divided by the number of matlab workers has to be an integer')
end

% Ricker wavelet (use f0=0 and t0=0 for an impulse)
model.f0 = 25;      %preak frequency
model.t0 = 0;      %time shift

%discretization options
%params.nppw    = 6;   %number of grid points per wavelength (should be [6 - 12])
%params.nwl_pml = 1; %physical size of the pml in terms of the longest wavelength in the start model

%computational options
params.nthreads = 1; %threads per spmd block (this is the lowest level of paralellism, used for the factorization, but not for forward/backward substitution)

% source/receiver grid and locations
%make sure these coordinates (in meters) are actually inside the
%computational domain!
model.zsrc = 25;
model.xsrc = 25:25:5725;
model.zrec = 200 ;
model.xrec = 25:25:5725;

model.nrec = length(model.xrec);
model.nsrc = length(model.xsrc);

% source matrix
Q = speye(length(model.xsrc));

%% startup
startup

current_dir = pwd;
addpath([current_dir,'/mbin']);
if exist([current_dir, resultsdir],'dir')==0
    mkdir([current_dir, resultsdir]);
end
cd results;

%% load and resample true velocity model in (m/s)
v = load(strcat(current_dir,'/data/cv.mat')); v=v.Data;

%use only a piece of the model
v = v(:,1535:1765);
v_true = v;
%----------
% Create background velocity model used as start model
v0_1=v(1:37,:);
v0_2=v(38:end,:);
[m2, n2]=size(v0_2);
S = opKron(opSmooth(size(v0_2,2),200),opSmooth(size(v0_2,1),200));
v0_2 = S*v0_2(:);
v0=[v0_1;reshape(v0_2,m2,n2)];

model.d = [6 25]; %z and x axis increment in the original model
model.o = [0 0];
model.n = [size(v,1) size(v,2)];
[model.z, model.x]=odn2grid(model.o,model.d,model.n);   %original grid

d_stable  = helm_stable_d( v,'m/s',PDEopts.HELM2D_CHEN9P,max(model.freq)+2);
if min(d_stable)<max(model.d)
    d_stable(3)=1;
    model.n(3)=1;
    model.d(3)=1;
    model.o(3)=0;
    min_v=min(v(:));
    max_v=max(v(:));
    min_v0=min(v0(:));
    max_v0=max(v0(:));
    
    [f2c,c2f,n_sub] = fine2coarse(model.n,model.d,d_stable,'cubic');
    v=f2c*v(:);
    v0=f2c*v0(:);
    
    v(v>max_v)=max_v;
    v(v<min_v)=min_v;
    
    v0(v0>max_v0)=max_v0;
    v0(v0<min_v0)=min_v0;
    
    model.d=d_stable;
    model.n=n_sub;
end

z_orig = model.z;
x_orig = model.x;

v=reshape(v,model.n);

% plot velocity model to be used
figure(1);
subplot(2,1,1);
imagesc(model.x,model.z,v);axis equal tight;colorbar;
xlabel('x [m]');ylabel('z [m]');title('True velocity model [m/s]');
subplot(2,1,2);
imagesc(model.x,model.z,reshape(v0,model.n));axis equal tight;colorbar;
xlabel('x [m]');ylabel('z [m]');title('start velocity model [m/s]');

model.n(3)=1;
model.d(3)=1;
model.o(3)=0;

% gridded slowness-squared values [km^2/s^2]
m = 1e6./v(:).^2;
m0= 1e6./v0(:).^2;

%% Set up the generation of synthetic data in the true velocity model
% load([current_dir datadir datafile]);
% D=Ds;clear Ds;
% D = reshape(D,model.nrec,model.nsrc,97); %73 is hardcoded, because it corresponds to the loaded data file. Change this number if you generate you own data
% D = D(:,:,1:model.nfreq);
% D=D(:);

%the lines below can be used to generate data and to add noise
 D           = F(m,Q,model,params);
 if save_data==1; Ds=gather(D); save('D_6ppw_5_29Hz','Ds'); clear Ds;  end;
 
% %add noise?
if add_noise==1;
    
    noise = randn(size(D,1),1)+1i*randn(size(D,1),1);
    noise = noise./norm(noise);
    noise = noise*NSR*norm(D(:));
    D     = D+noise;

if save_data==1; Ds=gather(D); save('D_noisy','Ds'); clear Ds;    end;
end

%% general setup
    if isfield(params,'pdefunopts')==0
        params.pdefunopts = PDEopts2D();
        params.pdefunopts.helm_pml=[];
    end
    if isfield(params,'lsopts')==0
        solve_opts = LinSolveOpts();
        solve_opts.solver = LinSolveOpts.SOLVE_LU;
        params.lsopts = solve_opts;
    end

%% Wavefield Reconstruction imaging

D=distributed(D);
params.wri       = true;     %use WRI. FWI is used when set to false
params.hessian   = 'sparse'; %type of Hessian approximation
params.lambda    = 1e2;      %this is the penalty parameter for WRI
fh = misfit_setup(m0,Q,D,model,params); %function handle which will return objective function, gradient and Hessian for the selected method (WRI/FWI)

[f_wri,g_wri] = fh(m0);
g_wri=gather(g_wri);

%% Reverse Time Migration (in frequency domain)

params.wri       = false;     % FWI (RTM) is used when set to false
params.extend    = false;
    
fh = misfit_setup(m0,Q,D,model,params); %function handle which will return objective function, gradient and Hessian for the selected method (WRI/FWI)
[f_fwi,g_fwi] = fh(m0);
g_fwi=gather(g_fwi);    

%% Plotting 
%(and some very basic image processing)
close all
n=model.n;
z=model.z;
x=model.x;

%plot true and background velocity
v0 = reshape(1e3./sqrt(m0),n); %convert back to velocity in [m/s]
figure(1)
subplot(2,1,1);set(gca,'Fontsize',14)
imagesc(x,z,v,[1500 4500]);colorbar;title('True model')
xlabel('x [m]');ylabel('z [m]');axis equal tight;h = colorbar;ylabel(h, 'Velocity [m/s]','FontSize',14);
subplot(2,1,2);set(gca,'Fontsize',14)
imagesc(x,z,v0,[1500 4500]);colorbar;title(['Initial model'])
xlabel('x [m]');ylabel('z [m]');axis equal tight;h = colorbar;ylabel(h, 'Velocity [m/s]','FontSize',14);

%plot migration results
gr = reshape((g_fwi),n);      % RTM
gp = reshape((g_wri),n);      % WRI method image

save('g_wri','g_wri');
save('g_fwi','g_fwi');

%optional: take the vertical diff to get rid of constant trends in the
%result
gr=(diff(gr));
gp=(diff(gp));

figure(2)
subplot(2,1,1);set(gca,'Fontsize',14)
imagesc(x,z,(gr),[-2e3 2e3]);colormap prx;title('Result RTM')
xlabel('x [m]');ylabel('z [m]');axis equal tight;
subplot(2,1,2);set(gca,'Fontsize',14)
imagesc(x,z,(gp),[-5e1 5e1]);colormap prx;title(['Result WRI, \lambda=',num2str(params.lambda)])
xlabel('x [m]');ylabel('z [m]');axis equal tight;

%save figures
if save_figures==1
    print(1,'-depsc',['True_and_inital_model']);
    print(2,'-depsc',['Migration_results']);
end

%% References
%[1]http://dx.doi.org/10.1093/gji/ggt258 Tristan van Leeuwen, Felix J. Herrmann, Geophysical Journal International,2013. Mitigating local minima in full-waveform inversion by expanding the search space.
%
%[2]https://slim.gatech.edu/content/penalty-method-pde-constrained-optimization Tristan van Leeuwen, Felix J. Herrmann. 2013. A penalty method for PDE-constrained optimization.
%
%[3]https://slim.gatech.edu/content/examples-penalty-method Bas Peters, Felix J. Herrmann, Tristan van Leeuwen. EAGE, 2014. Wave-equation based inversion with the penalty method: adjoint-state versus wavefield-reconstruction inversion.

 
