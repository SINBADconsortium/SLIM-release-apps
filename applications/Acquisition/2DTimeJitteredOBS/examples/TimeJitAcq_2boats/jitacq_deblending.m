% This script simulates the time-jittered (blended) marine acquisition scenario, and recovers conventional 
% data by sparse inversion (via one-norm minimization).
%
% Acquisition scenario: time-jittered OBC marine acquisition with two source vessels and two airgun arrays each
% (i.e., airgun arrays fire at jittered instances in time which translate to jittered shot locations for 
% a given speed of the source vessel)
%
% The results are stored in the path defined in the script setpath.m


% Set the parameters (by running the corresponding params.m script)
% NOTE: all the parameters are saved in a .mat file in the corresponding expdir
jitacq_deblending_params;


% Load seismic data
D = ReadSuFast(fname_input);
D = reshape(D, nt_org, nr_org, ns_org);
D = D(1:nt, 1:nr, 1:ns);


% Setup the time-jittered (or blended) acquisition
jitacq = jitter_airgunarrays(ns, ds, dt, p, nboats, randseed, boatspeed, tfireint_min, tdelay, delayboat, fig);

% Dimensions for time-jittered acquisition
jitdim = [(length(jitacq.tfirejitgrid) + nt - 1)  nr];

% Sampling operator 
if nboats == 1
   RM = opJitTimeShot1boat2arrays([nt nr ns], jitdim, jitacq);
elseif nboats == 2
   RM = opJitTimeShot2boats2arrays([nt nr ns], jitdim, jitacq);
end

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


% Setup the measurement operator [A], and generate measurements [b]
% 2D curvelets are applied on the receiver-source plane
% opCurvelet(M, N, nbscales, nbangles, finest, ttype, is_real)
C = opCurvelet(nr, ns, max(1,ceil(log2(min(nr,ns)) - 3)), 16, 1, 'ME', 0);

% 1D (spline) wavelets are applied on the time-axis 
% opSplineWaveletSPOT(M, N, filt, smooth, levels)
% NOTE: nt must be a power of 2 (required by the FFT operation performed within this operator) 
W = opSplineWaveletSPOT(nt, 1, nt, 3, 5);                   

% 2D curvelets Kroneckered with 1D wavelets act as a sparsifying transform
% oppKron2Lo: Kronecker tensor product to act on a distributed vector
S = oppKron2Lo(C, W', 1);                                      

% Measurement operator
A = RM*S'; 
 
% Measurements
b = jitD(:);


% Solve the one-norm minimization problem
% spgl1(A, b, tau, sigma, x, options)
% xest: estimated (synthesis) curvelet coefficients
% ctime: computation time
tic; xest = spgl1(A, b, spgl1_tau, spgl1_sigma, spgl1_x, options); ctime = toc;


% Recover deblended data [Dest], compute residual [Ddiff] and signal-to-noise ratio [SNR]
Dest = reshape(S'*xest, [nt nr ns]);
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

