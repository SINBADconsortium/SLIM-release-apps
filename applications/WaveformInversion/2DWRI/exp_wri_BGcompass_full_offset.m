%% Wavefield Reconstruction Inversion
%this is a synthetic example of WRI using the true source wavelet. Sources
%are close to the water surface and receivers are at the water bottom. Full
%source and reveir offset is used.

% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: July, 2014. Last updated February 2016.

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)

%This script will generate data and invert using the same modeling code.
%Model is the Compass model (part of it). Correct source function
%is used. 



%% startup
close all
startup

%data folder, files and results directory
basedir    = pwd;
resultsdir = [basedir,'/results/'];
datadir    = [basedir,'/data/'];
%datafile   = 'BGcompass_obn_fulloffset_5_29Hz';
if ~isdir(resultsdir)
    mkdir(resultsdir);
end
cd(resultsdir)

%% Main preferences (all user input is defined in this section)
label              = 'workspaceWRI_BG_5_9HZ_2sweep'; 
label_intermediate = 'WRI_BG_5_9_HZ_2sweep';

params.coarsen_grid = 0;

inversion_algorithm = 'WRI';
%inversion_algorithm = 'FWI';

params.lambda    = 1e2;      %this is the penalty parameter for WRI

%choose optimization algorithm: (see /SLIM-release-apps/applications/WaveformInversion/2DWRI/mbin/multi_stage_waveform_inversion.m for examples)
optimization_algorithm = 'TMP_lbfgs';        % Two-Metric-Projecton lbfgs

% solver: direct or iterative
iterative_solver = 0;
direct_solver    = 1;

%if iterative: use the lines below, if commented out, a direct solver will
%be used
if direct_solver==1 && iterative_solver ==1
    error('error:select iterative solver or direct solver, not both')
elseif iterative_solver ==1
    params.pde_solve = 'iter';
    params.ls_solver = 'solve_cgmn';
    params.precond   = 'prec_carp';
    params.refine    = 0;
    params.ls_iter   = inf;
    params.ls_tol    = 1e-4;
    params.sim_rec   = 0;
    params.n_sim_rec = 10;
elseif direct_solver ==1
    %params.nthreads = 4; %threads per spmd block (this is the lowest level of paralellism, used for the factorization)
    if isfield(params,'pdefunopts')==0
        params.pdefunopts = PDEopts2D();
        params.pdefunopts.helm_pml=[];
    end
    if isfield(params,'lsopts')==0
        solve_opts = LinSolveOpts();
        solve_opts.solver = LinSolveOpts.SOLVE_LU;
        params.lsopts = solve_opts;
    end
end

%logistics
save_vars.data           = 0;
save_vars.workspace      = 1;
save_vars.model_iters    = 1; %save model estimate at every iteration
save_vars.model_batch    = 1; %save the model estimate at the end of every frequency batch


% Gaussian noise
add_noise    = 0;
NSR          = 0.1; %noise to signal ration (2norm) ; NSR = norm(noise)/norm(signal)

%use simultaneous sources (yes=1). Redrawing means new random source
%weights are selected at the end of every iteration. Redrawing is mostly suitable for
%gradient-descent or Gauss-Newton type algorithms, but not for quasi-Newton
%algorithms.
params.sim_src         = 0;
params.n_sim_src       = 20;
params.redraw_sources  = 0;
params.redraw_interval = inf; %redraw every xx iterations of the optimization algorithm
%note: simultaneous sources without redrawing are still redrawn at the
%beginning of each frequency batch
params.get_new_fgH_redraw = 0; %obtains new f g and H after redrawing for reduced space methods (this takes extra PDE solves, but often makes the algorithm converge faster overall and to lower function values)

params.C        = 1; %this is a source/receiver mask in case of marine-aquisition, otherwise it is 1 (for land aquisition or ocean-bottom node)
params.srcest   = 0; %use source function estimation

%define frequency batches to be used (one at a time)
nbatch=4;
for i=1:nbatch
    f_a{i} = linspace(5+(i-1),6+(i-1),4);
end
f_a{end}=f_a{end}(1:3);
%pass through the data again ( 5->6->...->28->29->5->6->...->28->29Hz)
f_a = [f_a f_a]; 

%or alternatively, go back in frequency:
%( 5->6->...->28->29->28->27->...->6->5Hz)
%f_a = [f_a fliplr(f_a)];

model.freq  = unique(cell2mat(f_a));
model.nfreq = length(model.freq); nfreq=length(model.freq);

spmd;nl=numlabs;end;
if mod(model.nfreq,nl{1})~=0
    error('number of frequencies divided by the number of matlab workers has to be an integer')
end

%set number of iterations (of the nonlinear optimization algorithm) in each frequency batch (can vary per batch)
params.it_per_batch(1)             = 20;
params.it_per_batch(2:length(f_a)) = 15; 

% Ricker wavelet (use f0=0 and t0=0 for an impulse)
model.f0 = 10;      %preak frequency
model.t0 = 0;    %time shift

% source/receiver grid and locations
%make sure these coordinates (in meters) are actually inside the
%computational domain!
model.zsrc = 25;
model.xsrc = 25:50:5725;
model.zrec = 200;
model.xrec = 25:50:5725;

%% Set up Constaints 

%bound constraint options (all optimization methods use bound constraints)
constraint.v_max   = 5000; %upper bound
constraint.v_min   = 1450; %lower bound

constraint.water_max=1525;
constraint.water_min=1475;

constraint.water_depth=204;

constraint.implementation = 'single_constraint';
constraint=constraint_defaults(constraint);
constraint.use_bounds       = 1;

%% load and resample true velocity model in (m/s)
v = load([datadir, 'cv.mat']); v=v.Data;

%use only a piece of the model
v = v(:,1535:1765);
v_true = v;

% Create background velocity model used as start model
v0_1=v(1:37,:);
v0_2=v(38:end,:);
[m2, n2]=size(v0_2);
S = opKron(opSmooth(size(v0_2,2),200),opSmooth(size(v0_2,1),2000));
v0_2 = S*v0_2(:);
v0=[v0_1;reshape(v0_2,m2,n2)];
v0=repmat(mean(v0,2),1,size(v0,2));


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

%plot velocity model to be used
figure;
subplot(2,1,1);
imagesc(model.x,model.z,v,[1500 4500]);axis equal tight;colorbar;colormap jet
xlabel('x [m]');ylabel('z [m]');title('True velocity model [m/s]');
subplot(2,1,2);
imagesc(model.x,model.z,reshape(v0,model.n),[1500 4500]);axis equal tight;colorbar;colormap jet
xlabel('x [m]');ylabel('z [m]');title('start velocity model [m/s]');

model.n(3)=1;
model.d(3)=1;
model.o(3)=0;

% gridded slowness-squared values [km^2/s^2]
m = 1e6./v(:).^2;
m0= 1e6./v0(:).^2;

%% Set up source and receiver arrays
model.nrec = length(model.xrec);
model.nsrc = length(model.xsrc);

% source matrix
Q = speye(length(model.xsrc));

%% Set up the generation of synthetic data in the true velocity model
%check if data file already exists, otherwise run script to generate the data
% if exist([datadir,datafile,'.mat'], 'file') == 2
%     load([datadir,datafile,'.mat']);
% else
%     warning('data file not found, starting data generation')
%     cd(basedir)
%     generate_data_full_offset_obn.m
%     load([datadir,datafile,'.mat']);
%     disp('succesfully generated and loaded data')
% end

% D=Ds;clear Ds;
% D = reshape(D,data_info.nrec,data_info.nsrc,length(data_info.freq));
% D=D(:);
D           = F(m,Q,model,params);

% %add noise?
if add_noise==1;
    noise = randn(size(D,1),1)+1i*randn(size(D,1),1);
    noise = noise./norm(noise);
    noise = noise*NSR*norm(D(:));
    D     = D+noise;
    if save_data==1; Ds=gather(D); save('D_noisy','Ds'); clear Ds;    end;
end

%% Run multi-stage WRI or FWI
[mk,T,f] = multi_stage_waveform_inversion(m0,D,Q,model,f_a,inversion_algorithm,optimization_algorithm,params,label,label_intermediate,save_vars,constraint);

if save_vars.workspace==1;
    clear D
    save(['workspace_',label],'-append','v_true','v0');
end;
