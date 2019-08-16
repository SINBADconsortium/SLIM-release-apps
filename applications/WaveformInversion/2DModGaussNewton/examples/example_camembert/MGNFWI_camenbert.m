% This example use a small model to test the FWI algorithum
% -----------------------------------------------
%
% In this example, we process FWI from low frequency band to high frequency
% band. 10 GN subproblems are solved for each frequency band which contains 
% 10 frequencies. For each GN subproblem, we use 7 randomly selected simultaneous
% shots and roughly 20 L1 solver iterations
% 
% We sugguest you run this script with parallel matlab, 10 workers will be the best
%
% Results will be saved as a mat file under results directory.
% ------------------------------------------------
%
% setup a random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
close all;

%% camenbert example

%% define velocity model

% grid
z = 0:10:1000;
x = 0:10:1000;
[o,d,n] = grid2odn(z,x);
[zz,xx] = ndgrid(z,x);

% background velocity 2500 [m/s]
m0 = 2500 + 0*xx;

% circular perturbation with radius 250 m and strenth 10%
dv = 0*xx; dv((xx-500).^2 + (zz-500).^2 <= 250^2) = .1*2500;

% plot
m  = m0 + dv;
figure;imagesc(x,z,m);xlabel('x [m]');ylabel('z [m]');title('velocity [m/s]');colorbar;
m0 = m0(:);
m  = m(:);

%% setup parameters for reflection experiment
model.minf        = 2.9; % minimal frequency (Hz)
model.maxf        = 22;  % maximal frequency (Hz)
model.sourcet     = 1; % source type: 1 sim source; 2 sequential source.
model.name        = ['/example_camenbert/FWIresult_camenbert']; % name of result
model.vmin  = 2000;model.vmax  = 3000;% min vel and max vel for the projection
model.nf          = 10; % number of frequencies in each frequency band
model.snf         = 10; % number of frequencies for each GN updates 
model.ol          = 5;  % number of overlap frequencies for two adjacent bands
model.gl          = 10; % grid length
model.xbound      = .1; % x dimension boundary size, percentage
model.zbound      = .1; % z dimension boundary size, percentage
model.nshot       = 25;% number of shot positions
model.nsim        = 3; % Number of sim-shots for each GN subproblem
model.sdep        = 30; % source depth, unit: meters
model.sp          = round(linspace(1,101,model.nshot)); % source position
model.rdep        = 20; % receiver depth, unit: meters
model.nrec        = 51; % number of receivers
model.rp          = round(linspace(1,101,51)); % receiver positions
model.water       = 3; % estimated water depth, unit: number of grid

% setup optimazation parameters
opts = Set_Paras( ...
'miter',1, ...		                  % iteration number for frequency band.
'mlinm',10, ...			 	  		  % iteration number for relinearization for each freq band
'sptrans','no', ...                   % process in the physical domain
'iterations',2, ... 				  % iteration number for L1 solver ...
'dispresult',1);   


%% simulate observation data of reflection data

% define a wavelet
model.t0 = .2;
model.f0 = 12;
model.freq = linspace(model.minf,model.maxf,30);

% Setup forward modeling parameters.(do NOT modifed this part)
model.o = [0 0 0];
model.d = [d 1];
model.n = [n 1];        
model.nb = [floor(model.xbound*n(2)) floor(model.zbound*n(1)) 0];
model.zsrc = model.sdep * model.gl; 
model.xsrc = (model.sp - 1) * model.gl ;
model.zrec = model.rdep * model.gl; 
model.xrec = (model.rp - 1) * model.gl ;
Q = (speye(model.nshot)./(model.gl^2));
Dobs  = gather(F(1e6./m.^2,Q,model)); 
Dobs  = reshape(Dobs,model.nrec,model.nshot,length(model.freq));

% wavelet information
wavelet.t0 = model.t0;
wavelet.f0 = model.f0;
wavelet.faxis = model.freq;




%% inversion

[results] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts);



%% setup modeling parameters for transmission experiment
model.rdep = 950;
model.name = ['/example_camenbert/FWIresult_camenbert_transmission']; % name o
model.zrec = model.rdep;
Dobs  = gather(F(1e6./m.^2,Q,model)); 
Dobs  = reshape(Dobs,model.nrec,model.nshot,length(model.freq));


%% inversion

[results] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts);




