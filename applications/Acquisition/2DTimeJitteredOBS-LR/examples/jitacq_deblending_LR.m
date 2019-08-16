% This script simulates the time-jittered (blended) marine acquisition scenario, and recovers conventional 
% data by rank-minimization.
%
% Acquisition scenario: time-jittered OBC marine acquisition with one source vessel and two airgun arrays
% (i.e., airgun arrays fire at jittered instances in time which translate to jittered shot locations for 
% a given speed of the source vessel)
%
% The results are stored in the path defined in the script setpath.m


% Set the parameters
jitacq_deblending_LR_params;

% Load seismic data
D = ReadSuFast(fname_input);
D = reshape(D, nt_org, nr_org, ns_org);
D = D(1:nt, 1:nr, 1:ns);

% Setup the time-jittered (or blended) acquisition
jitacq = jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig);

% Dimensions for time-jittered acquisition
jitdim = [(length(jitacq.tfirejitgrid) + nt - 1)  nr];

% Frequency samples
wn = [0:nt/2,-nt/2+1:-1]*2*pi/(nt*dt);

% Sampling operator (Actual sampling operator)
RM = opJittshift([nt nr ns], jitdim, jitacq, wn');

% Test the operator for the dottest
if test, dottest(RM, 1, 'verbose'); end

% Generate time-jittered (or blended) data volume [jitD]
jitD = reshape(RM*D(:), jitdim);

% Save jittered data volume
write_ACQdata(fname_jitdata, jitD, dt, dr, ds)


% Recover deblended data by applying the adjoint of the sampling operator
Dest_adj = reshape((RM')*jitD(:), [nt nr ns]);
 
% Save (adjoint) recovered data
write_ACQdata(fname_adjrecov, Dest_adj, dt, dr, ds)
 
% Measurements
b = jitD(:);

% spgl1 parameter
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                   'decTol',1e-6,...
                   'project', @TraceNorm_project_TJM, ...
                   'primal_norm', @TraceNorm_primal_TJM, ...
                   'dual_norm', @TraceNorm_dual_TJM, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations',iteration,...
                   'fid',fid,...
                   'weights', []);
params.nm   = nr;
params.nh   = 2*ns-1;
params.nf   = floor(nt/2)+1;
params.nr   = nr;
params.nc   = ns;
params.nt   = nt;
params.mode = 1;
params.RM   = RM;
params.funForward = @Time_Jitt_LR;
params.ls   = 1;
params.mode = 1;
sigma       = norm(jitD,2);
params.k    = rank;
LInit       = randn(params.nf,params.nm,params.k)+1i*randn(params.nf,params.nm,params.k);
RInit       = randn(params.nf,params.nh,params.k)+1i*randn(params.nf,params.nh,params.k);
xinit       = initiweight*[vec(LInit);vec(RInit)];
tau         = norm(xinit,1);
sigmafact   = sigmaerror*sigma;

% Solve the rank minimization problem
tic; xest = spgl1(@Time_Jitt_LR,b,tau,sigmafact,xinit,opts,params); ctime = toc;

% Recover deblended data [Dest], compute residual [Ddiff] and signal-to-noise ratio [SNR]
Ft   = opDFTR(params.nt);
Ft   = opKron(opDirac(params.nr*params.nc),Ft);
L    = xest(1:params.nf*params.nm*params.k);
R    = xest(params.nf*params.nm*params.k+1:end);
L    = reshape(L,params.nf,params.nm,params.k);
R    = reshape(R,params.nf,params.nh,params.k);
MH   = opMH(params.nr,params.nc);
Dest = zeros(params.nf,params.nr*params.nc);

% tranform tehdata from MH domain to SR domain
for i = 1:params.nf
    Dest(i,:) = MH'*vec(squeeze(L(i,:,:))*squeeze(R(i,:,:))');
end

% perform the inverse-FFT along the spatial axis
Dest = Ft'*vec(Dest);
Dest = reshape(Dest,nt,nr,ns);
Ddiff = D - Dest;
SNR = -20*log10(norm(Ddiff(:))/norm(D(:)));

% Save recovered data
write_ACQdata(fname_L1recov, Dest, dt, dr, ds)

% Save difference data
write_ACQdata(fname_L1diff, Ddiff, dt, dr, ds)

% Save ctime and SNR to .mat file
save('comptime_SNR.mat', 'ctime', 'SNR');


% Return to the examples directory
cd(curdir)

