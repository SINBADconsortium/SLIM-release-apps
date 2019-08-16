% This script simulates the over/under, blended marine data and deblends the data into its constituent source
% components by promoting low rank via nuclear-norm minimization.
%
% Acquisition scenario: over/under (blended) marine acquisition with one source vessel and two airgun arrays
% (one airgun array at a depth of 10.0 m and the other at 15.0 m)
%
% The results are stored in the path defined in the script setpath.m

% Open matlabpool
% NOTE: use this when running the code interactively
parpool(4);

% Set the parameters 
% NOTE: all the parameters are saved in a .mat file in the corresponding expdir
SourceSep_params;

% Load seismic data
D1 = reshape(ReadSuFast(fname_input1,'b'), nt, nr, ns);
D2 = reshape(ReadSuFast(fname_input2,'b'), nt, nr, ns);

% Generate blended data
[delay_s2, D2_shift, D_blend] = gen_blended_data(rseed, dt, nt, ns, nr, D1, D2);

% Save shifted data (source 2)
save_data('data_TRS', fname_D2_shift, D2_shift, dt)

% Save blended data (in time domain)
save_data('data_TRS', fname_D_blend_TRS, D_blend, dt)

% Transform data into frequency domain
D1_FRS = tfdomain(D1, nt, nfreq, nr, ns, 1);
D2_FRS = tfdomain(D2, nt, nfreq, nr, ns, 1);
D2_shift_FRS = tfdomain(D2_shift, nt, nfreq, nr, ns, 1);
D_blend_FRS = tfdomain(D_blend, nt, nfreq, nr, ns, 1);

% Save blended data (in frequency domain)
save_data('data_FRS', fname_D_blend_FRS, permute(D_blend_FRS,[3 1 2]), dfreq)

% Time delay matrix
params.delay = delay_s2;
time = params.depth/params.watervel;
tdelay = time + repmat(params.delay,ns,1);

% Range of frequencies with significant energy
ind = round(freq_min/dfreq)+1 : round(freq_max/dfreq)+1;
freq_range = freq_axis(ind);

% Set value of rank for each frequency slice
rk = floor(linspace(rk_min, rk_max, length(freq_range)));

% HSS level indices
HSS_index = HSS_cell(ones(nr,ns), HSS_level);
HSS_index = distributed(HSS_index);

% MAIN EXECUTION LOOP

% TRS ---> Time-Receiver-Source
% FRS ---> Frequency-Receiver-Source
% SR (or sr) ---> Source-Receiver domain
% MH (or mh) ---> Midpoint-Offset domain
% S1est ---> recovered/separated frequency slice for source 1 in the mh domain 
% S2est ---> recovered/separated frequency slice for source 2 in the mh domain 
% S1estsr ---> recovered/separated frequency slice for source 1 in the sr domain 
% S2estsr ---> recovered/separated frequency slice for source 2 in the sr domain 
% D1recov ---> recovered/separated data cube for source 1 in the FRS domain 
% D2_shiftrecov ---> recovered/separated data cube for source 2 in the FRS domain 

% Number of columns 
nc = ns;

% Initialize recovered data cubes
D1_recov = zeros(nfreq,nr,nc);
D2_shiftrecov = zeros(nfreq,nr,nc);

for j = 1:length(freq_range)

    % Initialize recovered matrices
    S1estsr = zeros(nr,nc);
    S2estsr = zeros(nr,nc);
    
    cur_rk = rk(j);
    omega = 2*pi*freq_range(j);
    b_input = squeeze(D_blend_FRS(j,:,:));

    spmd
        
        HSS_index_local = getLocalPart(HSS_index);
        
        for k = 1:size(HSS_index_local,2)

	    params.rk = cur_rk;
            params.omega = omega;

            % Perform HSS partitioning on the input frequency slice (in SR domain)            
            b = HSS_marine(b_input,nr,nc,HSS_index_local(:,k),1);
            [nrsub,ncsub] = size(b);
            params.sign = ones(size(b));

            % Midpoint-offset operator  
            MH = opMH(nrsub,ncsub);
            params.mhnumr = min(nrsub,ncsub) + floor(abs(nrsub-ncsub)/2);
            params.mhnumc = nrsub + ncsub - 1;

            % Transform each HSS block from SR to MH domain     
            b = reshape(MH*vec(b),params.mhnumr,params.mhnumc);
            params.sign = reshape(MH*vec(params.sign),params.mhnumr,params.mhnumc);
            [nmsub,nhsub] = size(b);
            params.tdelay_sub = HSS_marine(tdelay,nr,nc,HSS_index_local(:,k),1);
            params.tdelay_sub = reshape(MH*vec(params.tdelay_sub),params.mhnumr,params.mhnumc);

            % Initial guess for L, R factors
            L1 = randn(params.mhnumr,params.rk) + 1i*randn(params.mhnumr,params.rk);
            R1 = randn(params.mhnumc,params.rk) + 1i*randn(params.mhnumc,params.rk);
            L2 = randn(params.mhnumr,params.rk) + 1i*randn(params.mhnumr,params.rk);
            R2 = randn(params.mhnumc,params.rk) + 1i*randn(params.mhnumc,params.rk);

            % SPGL1 
            % NOTE: input b is in the MH domain
            tau = norm(xinit,1);
            sigma = norm(b(:),2);
            sigmafact = 1e-3*sigma;
            xinit = 1e-6*[vec(L1);vec(R1);vec(L2);vec(R2)];
            [output,r,g,info] = spgl1JULY(@NLfunForward_parallel, b(:), tau, sigmafact, xinit, opts, params);

            % Compute recovered L, R factors
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

            % Compute separated frequency slices (in MH domain) for each HSS block
            S1est = L1*R1';
            S2est = L2*R2';
          
            % MH to SR domain for each HSS block
            S1estsr_sub = reshape(MH'*vec(S1est),nrsub,ncsub);
            S2estsr_sub = reshape(MH'*vec(S2est),nrsub,ncsub);

            % Concatenate all HSS blocks            
            S1estsr = S1estsr + HSS_marine(S1estsr_sub,nr,nc,HSS_index_local(:,k),-1);
            params.tdelay_sub = reshape(MH'*vec(params.tdelay_sub),nrsub,ncsub);
            S2estsr_sub = exp(-1i*params.omega*(params.tdelay_sub)).*S2estsr_sub;
            S2estsr = S2estsr + HSS_marine(S2estsr_sub,nr,nc,HSS_index_local(:,k),-1);

        end
       
        % Global summation (for spmd only)  
        S1estsr = pSPOT.utils.global_sum(S1estsr);
        S2estsr = pSPOT.utils.global_sum(S2estsr);

    end

    % Store the recovered (or separated) frequency slices 
    D1_recov(j,:,:) = S1estsr{1};
    D2_shiftrecov(j,:,:) = S2estsr{1};
    
    % Compute SNRs 
    SNR1(j) = -20*log10(norm(squeeze(D1_FRS(j,:,:)) - S1estsr{1},'fro')/norm(squeeze(D1_FRS(j,:,:)),'fro'));
    SNR2(j) = -20*log10(norm(squeeze(D2_shift_FRS(j,:,:)) - S2estsr{1},'fro')/norm(squeeze(D2_shift_FRS(j,:,:)),'fro'));
    
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

% Close matlabpool 
% NOTE: use this when running the code interactively
delete(gcp('nocreate'));

% Return to the examples directory
cd(curdir)

