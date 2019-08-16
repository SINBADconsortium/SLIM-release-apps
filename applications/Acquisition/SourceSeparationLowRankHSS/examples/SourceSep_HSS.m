% This script simulates the over/under, blended marine data and deblends the data into its constituent source
% components by promoting low rank via nuclear-norm minimization.
%
% Acquisition scenario: over/under (blended) marine acquisition with one source vessel and two airgun arrays
% (one airgun array at a depth of 10.0 m and the other at 15.0 m)
%
% The results are stored in the path defined in the script setpath.m


% Set the parameters 
% NOTE: all the parameters are saved in a .mat file in the corresponding expdir
SourceSep_params;

% Load seismic data
D1 = reshape(ReadSuFast(fname_input1,'b'), nt, nr, ns);
D2 = reshape(ReadSuFast(fname_input2,'b'), nt, nr, ns);

% Generate blended data
[delay_s2, D2_shift, D_blend] = gen_blended_data(rseed, dt, nt, ns, nr, D1, D2);

% Save data for source 1
save_data('data_TRS', fname_D1, D1, dt)

% Save shifted data (source 2)
save_data('data_TRS', fname_D2_shift, D2_shift, dt)

% Save blended data (in time domain)
save_data('data_TRS', fname_D_blend_TRS, D_blend, dt)


% Transform data into frequency domain
t = 0:dt:(nt-1)*dt;
D1spec = tfdomain(D1, nt, [], nr, ns, 1);
D2spec = tfdomain(D2, nt, [], nr, ns, 1);
D2_shift_spec = tfdomain(D2_shift, nt, [], nr, ns, 1);
[D_blend_spec, dim] = tfdomain(D_blend, nt, [], nr, ns, 1);
freq_axis = dfreq*(0:dim(1)-1);
nfreq = dim(1);

% Save blended data (in frequency domain)
save_data('data_FRS', fname_D_blend_FRS, permute(D_blend_spec,[3 1 2]), dfreq)

% Time delay matrix
params.delay = delay_s2;
time = 2*params.depth/params.watervel;
tdelay = repmat(params.delay, ns, 1);

% Range of frequencies with significant energy
freq_range = find(freq_axis <= f_high);

% Set value of rank for each frequency slice
rk = floor(linspace(rk_min, rk_max, length(freq_range)));

% HSS level indices
HSS_index = HSS_cell(ones(nr,ns), level);


% MAIN EXECUTION LOOP

% TRS ---> Time-Receiver-Source
% FRS ---> Frequency-Receiver-Source
% SR (or sr) ---> Source-Receiver domain
% MH (or mh) ---> Midpoint-Offset domain
% S1est ---> estimated/recovered/separated frequency slice for source 1 in the mh domain 
% S2est ---> estimated/recovered/separated frequency slice for source 2 in the mh domain 
% S1estsr ---> estimated/recovered/separated frequency slice for source 1 in the sr domain 
% S2estsr ---> estimated/recovered/separated frequency slice for source 2 in the sr domain 
% D1recov ---> recovered/separated data cube for source 1 in the FRS domain 
% D2_shiftrecov ---> recovered/separated data cube for source 2 in the FRS domain 

% Initialize recovered data cubes
D1recov = zeros(nr,ns,length(freq_axis));
D2_shiftrecov = zeros(nr,ns,length(freq_axis));

for j = 1:length(rk)

   % Initialize recovered matrices (for source 1 (S1) and source 2 (S2) in sr domain)
   S1estsr = zeros(nr,ns);
   S2estsr = zeros(nr,ns);

   params.rk = rk(j);
   params.omega = 2*pi*freq_axis(j);
   b_input = squeeze(D_blend_spec(:,:,j));

   % LOOP OVER HSS BLOCKS
   for k = 1:size(HSS_index,2)
 
      b = HSS_macq(b_input,nr,ns,HSS_index(:,k),1);
      [nrsub,ncsub] = size(b);
      params.tdelay_sub = HSS_macq(tdelay,nr,ns,HSS_index(:,k),1);

      % MH operator
      MH = opMH(nrsub,ncsub);
      params.srnumr = nrsub;
      params.srnumc = ncsub;
      params.mhnumr = min(nrsub,ncsub) + floor(abs(nrsub-ncsub)/2);
      params.mhnumc = nrsub + ncsub - 1;   

      % RHS
      params.afun = @(x)afunacq_HSS(x,MH,params,1);
      params.afunT = @(x)afunacq_HSS(x,MH,params,-1);

      % Initial guess
      L1 = randn(params.mhnumr,params.rk) + 1i*randn(params.mhnumr,params.rk);
      R1 = randn(params.mhnumc,params.rk) + 1i*randn(params.mhnumc,params.rk);
      L2 = randn(params.mhnumr,params.rk) + 1i*randn(params.mhnumr,params.rk);
      R2 = randn(params.mhnumc,params.rk) + 1i*randn(params.mhnumc,params.rk);
      xinit = 1e-6*[vec(L1);vec(R1);vec(L2);vec(R2)];

      % SPGL1 
      tau = norm(xinit,1);
      sigma = norm(b(:),2);
      sigmafact = 1e-3*sigma;
      [output,r,g,info] = spgl1(@NLfunForward, b(:), tau, sigmafact, xinit, opts, params);

      % Compute L, R factors
      e1 = params.mhnumr*params.rk;
      e2 = params.mhnumc*params.rk;

      L1 = output(1:e1);
      R1 = output(e1+1:e1+e2);
      L1 = reshape(L1,params.mhnumr,params.rk);
      R1 = reshape(R1,params.mhnumc,params.rk);

      L2 = output(e1+e2+1:2*e1+e2);
      R2 = output(2*e1+e2+1:end);
      L2 = reshape(L2,params.mhnumr,params.rk);
      R2 = reshape(R2,params.mhnumc,params.rk);

      % Compute the separated frequency slices (in mh domain) for each HSS block
      S1est = L1*R1';
      S2est = L2*R2';

      % MH to SR domain (for each HSS block)
      S1estsr_sub = reshape(MH'*vec(S1est),nrsub,ncsub);
      S2estsr_sub = reshape(MH'*vec(S2est),nrsub,ncsub);

      % Concatenate all HSS blocks 
      S1estsr = S1estsr + HSS_macq(S1estsr_sub,nr,ns,HSS_index(:,k),-1);
      S2estsr_sub = exp(-1i*params.omega*(time + params.tdelay_sub)).*S2estsr_sub;
      S2estsr = S2estsr + HSS_macq(S2estsr_sub,nr,ns,HSS_index(:,k),-1);

   end

   % Store the recovered (or separated) frequency slices 
   D1recov(:,:,j) = S1estsr;
   D2_shiftrecov(:,:,j) = S2estsr;

end

% END OF MAIN LOOP


% Convert recovered data cubes from FRS to TRS domain
D1recov_TRS = tfdomain(D1recov, nt, nfreq, nr, ns, -1);
D2_shiftrecov_TRS = tfdomain(D2_shiftrecov, nt, nfreq, nr, ns, -1);

% Save data
save_data('data_FRS', fname_D1recov_FRS, permute(D1recov,[3 1 2]), dfreq)
save_data('data_TRS', fname_D1recov_TRS, D1recov_TRS, dt)
save_data('data_FRS', fname_D2recov_FRS, permute(D2_shiftrecov,[3 1 2]), dfreq)
save_data('data_TRS', fname_D2recov_TRS, D2_shiftrecov_TRS, dt)


% Return to the examples directory
cd(curdir)

