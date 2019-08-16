% Illustrates the benifit of solving constrained formulations of FWI.
% Various examples of constraints for a simple 2D model in a well-to-well
% setting. Receivers in one well, sources in the other. This code shows how
% to use the constrained optimization toolbox 
% (https://github.com/SINBADconsortium/SLIM-release-apps/tree/master/applications/WaveformInversion/ConstrainedFWI) 
% via a few simple examples.
%
% It should take a couple of minutes to run. To run
% in parallel, use parpool(6), because there are 6 frequencies.

% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Last update: September 2016

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)

close all
assert(parpool_size==6,'size of parallel pool should be 6 workers');

%% noise level (can be NSR=0, NSR : ||noise||_2 / ||signal||_2
NSR=0.5;

%% startup
disp('Do not forget to run starup.m located in this folder') 
current_dir = pwd;

%% Define general options for constraints
constraint.projection_accuracy = 'accurate';
[constraint]=constraint_defaults(constraint); %this gets the default settings that are required for the toolbox to work

%% optimization options, for the nonlinear solver (spectral) projected-gradient
options.maxIter     = 20;
options.memory      = 5;
options.testOpt     = 0;
options.progTol     = 10*eps;
options.fprogTol    = 10*eps;
options.useSpectral = 1;

%% define velocity model

% grid
z = 0:10:1000;
x = 0:10:1000;
[o,d,n] = grid2odn(z,x);
[zz,xx] = ndgrid(z,x);

% background velocity 2500 [m/s]
v0 = ones(n).*2500;

% square perturbation
v = v0;
center_x = round(length(x)/2); center_z = round(length(z)/2);
v(center_z-15:center_z+15,center_x-15:center_x+15)=2500+.1*2500;

vpmin=min(v(:))-100;
vpmax=max(v(:))+100;

% plot
figure;imagesc(x,z,v,[vpmin vpmax]);xlabel('x [m]');ylabel('z [m]');title('true velocity [m/s]');colorbar;colormap jet

%% set up the experiment: sources, receivers frequencies, source function

% grid
model_t.o = o; 
model_t.d = d;
model_t.n = n;

% frequencies in Hz.
model_t.freq = linspace(5,15,6); nfreq = length(model_t.freq);

% Ricker wavelet with peak frequency f0 and phase shift t0
model_t.f0 = 15;
model_t.t0 = 0;

model_t.xsrc = 0:100:1000; nsrc = length(model_t.xsrc);
model_t.zsrc = 50;
model_t.xrec = 0:10:1000;  nrec = length(model_t.xrec);
model_t.zrec = 950;
model_t.n(3)=1;
model_t.d(3)=1;
model_t.o(3)=0;

% source matrix, each column is a source function defined on the source
% grid [model.zsrc, model.xsrc].
Q = speye(nsrc);

% convert to squared slowness [km^2/s^2]
m = 1e6./v(:).^2;

% initial model
m0 = 1e6./v0(:).^2;

%some linear algebra options
params=[];
if isfield(params,'pdefunopts')==0
    params.pdefunopts = PDEopts2D();
end
if isfield(params,'lsopts')==0
    solve_opts = LinSolveOpts();
    solve_opts.solver = LinSolveOpts.SOLVE_LU;
    params.lsopts = solve_opts;
end

%% create synthetic 'observed' data
D_t = F(m,Q,model_t,params);

%add noise?
    noise = randn(size(D_t,1),1)+1i*randn(size(D_t,1),1);
    noise = noise./norm(noise);
    noise = noise*NSR*norm(D_t(:));
    D_t     = D_t+noise;


%% create function handle for misfit
%this is the only component that needs to be replaced if you want to use
%your own codes which provides a function value and gradient. It takes the
%velocity model (as a vector) as input, and ouputs the scalar function
%value (f) and gradient (g, as a vector)
%[f,g]=fh(velocity_model)

% params.wri = false;
% params.dist_mode = 'freq';
fh = misfit_setup(m0,Q,D_t,model_t,params);


%% use SPG with only bound constraints

constraint.use_bounds       = 1;
constraint.v_max   = vpmax; %upper bound
constraint.v_min   = vpmin; %lower bound

%obtain projectors
params.freq_index=[1:length(model_t.freq)];
P = setup_constraints_freq_FWI_2D(constraint,model_t,model_t,params,m0,1); %get projectors onto each set separately
funProj =@(input) P{1}(input); %only one constraint in this case
[mn_t,fsave,funEvals,projects,iter_save] = minConf_SPG(fh,m0,funProj,options);
%[mn_t,fsave,funEvals,projects] = minConf_PQN(fh,m0,funProj,options);
mn_t=1e3./sqrt(mn_t); %convert back to m/s

% plot
figure;imagesc(x,z,reshape(real(mn_t),n),[vpmin vpmax]);xlabel('x [m]');ylabel('z [m]');title('reconstruction from transmission data - bounds only');
colorbar;axis equal tight;colormap jet;


%% solve transmission experiment - bounds and Fourier based minimum smoothness

% define constraints here
[constraint]=constraint_defaults(constraint);
constraint.use_bounds       = 1;
constraint.use_min_smooth   = 1;

%bound constraint options (all optimization methods use bound constraints)
constraint.v_max   = vpmax; %upper bound
constraint.v_min   = vpmin; %lower bound

% minimum smoothness    full none aspect
constraint.smoothpars    = [3,  3,   1];

%obtain projectors
params.freq_index=[1:length(model_t.freq)];
P = setup_constraints_freq_FWI_2D(constraint,model_t,model_t,params,m0,1); %get projectors
funProj = @(input) Dykstra_cyclic(input,P,constraint.options_dyk); % get projector onto the intersection

% use SPG with constraints
[mn_t,fsave,funEvals,projects,iter_save] = minConf_SPG(fh,m0,funProj,options);
mn_t=1e3./sqrt(mn_t);

% plot
figure;imagesc(x,z,reshape(mn_t,n),[vpmin vpmax]);xlabel('x [m]');ylabel('z [m]');title('reconstruction from transmission data - bounds and Fourier minimum smoothness');
axis equal tight;colormap jet;colorbar

%% transmission experiment - Bound constraints and rank

% define constraints here
[constraint]=constraint_defaults(constraint);
constraint.use_bounds       = 1;
constraint.use_rank         = 1;

%bound constraint options (all optimization methods use bound constraints)
constraint.v_max   = vpmax; %upper bound
constraint.v_min   = vpmin; %lower bound

%rank
constraint.r_initial   = 2; %initial rank of the estimated model
constraint.r_increment = 0; %increase the rank of the estimated model by this number each new frequency batch

%obtain projectors
params.freq_index=[1:length(model_t.freq)];
P = setup_constraints_freq_FWI_2D(constraint,model_t,model_t,params,m0,1); %get projectors
funProj = @(input) Dykstra_cyclic(input,P,constraint.options_dyk); % get projector onto the intersection

% use SPG with constraints
[mn_t,fsave,funEvals,projects,iter_save] = minConf_SPG(fh,m0,funProj,options);
mn_t=1e3./sqrt(mn_t);

% plot
figure;imagesc(x,z,reshape(mn_t,n),[vpmin vpmax]);xlabel('x [m]');ylabel('z [m]');title('reconstruction from transmission data - bounds and rank constraints');
axis equal tight;colormap jet;colorbar

%% transmission experiment - Bound constraints and total-variation

[constraint]=constraint_defaults(constraint);
constraint.use_bounds       = 1;
constraint.use_TV           = 1;
constraint.v_max   = vpmax; %upper bound
constraint.v_min   = vpmin; %lower bound

%Total Variation
constraint.TV_ini_abs = 0.3443; %this is the true TV, to see the effect of under/over estimating the TV, just change this number
%constraint.TV_initial    = 1;
%constraint.TV_increment  = 1;
constraint.TV_units      = 's^2/m^2';

%obtain projectors
params.freq_index=[1:length(model_t.freq)];
P = setup_constraints_freq_FWI_2D(constraint,model_t,model_t,params,m0,1); %get projectors
funProj = @(input) Dykstra_cyclic(input,P,constraint.options_dyk); % get projector onto the intersection

% use SPG with constraints
[mn_t,fsave,funEvals,projects,iter_save] = minConf_SPG(fh,m0,funProj,options);
mn_t=1e3./sqrt(mn_t);

% plot
figure;imagesc(x,z,reshape(mn_t,n),[vpmin vpmax]);xlabel('x [m]');ylabel('z [m]');title('reconstruction from transmission data - bounds and TV constraints');
axis equal tight;colormap jet;colorbar

%% transmission experiment - Bound constraints and total-variation and rank

% define constraints here
[constraint]=constraint_defaults(constraint);
constraint.use_bounds       = 1;
constraint.use_rank         = 1;
constraint.use_TV           = 1;

%bound constraint options (all optimization methods use bound constraints)
constraint.v_max   = vpmax; %upper bound
constraint.v_min   = vpmin; %lower bound

%Total Variation
constraint.TV_ini_abs = 0.3443;
%constraint.TV_increment  = 1;
constraint.TV_units      = 's^2/m^2';

%rank
constraint.r_initial   = 2; %initial rank of the estimated model
constraint.r_increment = 0; %increase the rank of the estimated model by this number each new frequency batch

%obtain projectors
params.freq_index=[1:length(model_t.freq)];
P = setup_constraints_freq_FWI_2D(constraint,model_t,model_t,params,m0,1); %get projectors
funProj = @(input) Dykstra_cyclic(input,P,constraint.options_dyk); % get projector onto the intersection

% use SPG with constraints
[mn_t,fsave,funEvals,projects,iter_save] = minConf_SPG(fh,m0,funProj,options);
mn_t=1e3./sqrt(mn_t);

% plot
figure;imagesc(x,z,reshape(mn_t,n),[vpmin vpmax]);xlabel('x [m]');ylabel('z [m]');title('reconstruction from transmission data - bounds, TV and rank constraints');
axis equal tight;colormap jet;colorbar
