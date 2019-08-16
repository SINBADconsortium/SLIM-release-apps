% This script simulates the time-jittered blended marine acquisition scenario, and recovers conventional (deblended) data by sparse inversion via one-norm minimization.
%
% Acquisition scenario: time-jittered OBC marine acquisition with one source vessel and two airgun arrays (i.e., airgun arrays fire at jittered instances in time which translate to jittered shot locations for a given speed of the source vessel).
%
% The results are stored in the path defined in the script setpath.m.


% Set the parameters
% NOTE: all the parameters are saved in a .mat file in the corresponding expdir
TimeJitAcq_params;

% Load data
D = ReadSuFast(fname_input);
D = reshape(D, nt_org, nr_org, ns_org);
D = D(1:nt, 1:nr, 1:ns);

% Setup the jittered (or blended) acquisition
jitacq = jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams);

% Save jittered acquisition parameters
save('jitacqparams.mat', 'jitacq')

% Load jittered acquisition parameters
%load([expdir 'jitacqparams.mat'])

% Dimensions for jittered acquisition
jitdim = [(length(jitacq.tfirejitgrid) + nt - 1)  nr];

% Frequency samples
wn = [0:nt/2,-nt/2+1:-1]*2*pi/(nt*dt); 

% Time grid indices for the sampling operator (when applied on NFFT generated data)
sjitpos = ([jitacq.sjitb1arr1 jitacq.sjitb1arr2]);
sjitpos_sort = sort([jitacq.sjitb1arr1 jitacq.sjitb1arr2]);
idx = zeros(1,length(sjitpos));
for k = 1:length(sjitpos); idx(k) = find(sjitpos_sort(k) == sjitpos); end

% Sampling operator 
RM = opJitNFDCT1boat([nt nr length(idx)], jitdim, jitacq, idx, wn', nr);

% Test the operator for the dottest
if test, dottest(RM, 1, 'verbose'); end

% Generated jittered, blended data
newD = genJitData(nboats, jitacq, D, ds);
jitD = reshape(RM*newD(:), jitdim);

% Save irregular data generated generated via NFFT (output from genJitData.m)
write_ACQdata(fname_irregdata, newD, dt, dr, ds)

% Save irregular data generated generated via NFFT (output from genJitData.m)
write_ACQdata(fname_jitdata, jitD, dt, dr, ds)

% Run the main algorithm
Dest = TimeJitAcq_Deblend(label, jitacq, nt, nr, ns, ds, RM, jitD, opt, options);

% Compute residual
Ddiff = D - Dest; 

% Save recovered data
write_ACQdata(fname_L1recov, Dest, dt, dr, ds)

% Save difference data
write_ACQdata(fname_L1diff, Ddiff, dt, dr, ds)


% Return to the examples directory
cd(curdir)

