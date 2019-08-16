% This script simulates the over/under, blended marine data and deblends the data into its constituent source
% components by sparsity promotion.
%
% Acquisition scenario: over/under (blended) marine acquisition with one source vessel and two airgun arrays
% (one airgun array at a depth of 8.0 m and the other at 12.0 m)
%
% The results are stored in the path defined in the script setpath.m

% Set the parameters 
% NOTE: all the parameters are saved in a .mat file in the corresponding expdir
SourceSep_params;

% Load delay times 
delay_S1 = rsf_read_all(fname_delay_S1);
delay_S2 = rsf_read_all(fname_delay_S2);

% Load time-shifted data
D1_shift = rsf_read_all(fname_D1_shift);
D2_shift = rsf_read_all(fname_D2_shift);
D1_shift_FRS = rsf_read_all(fname_D1_shift_FRS);
D2_shift_FRS = rsf_read_all(fname_D2_shift_FRS);

% Load blended data
% FRS: data in frequency-receiver-source dimensions
D_blend = rsf_read_all(fname_D_blend);
D_blend_FRS = rsf_read_all(fname_D_blend_FRS);

% Use only one side of the FX spectra, since the other half is symmetric
if mod(nt,2) == 0 
  nfreq = floor((nt/2)+1);
else
  nfreq = ceil(nt/2);
end

% Curvelet operator
C = opCurvelet(nr, ns, max(1,ceil(log2(min(nr,ns)) - 3)), 16, 1, 'WRAP', 0);

% Set other parameters
params.nt = nt;
params.nr = nr;
params.ns = ns;
params.curvlength = size(C,1);
params.delay_S1 = delay_S1(231:end);
params.delay_S2 = delay_S2(231:end);

% Time delay matrix
params.tdelay_S1 = repmat(params.delay_S1,ns,1);
params.tdelay_S2 = repmat(params.delay_S2,ns,1);

% Indices of frequencies to separate
% NOTE - all frequencies used for separation 
ind = 1 : round(freq_max/dfreq)+1;

% NOTE: this code can be run in parallel by splitting the frequencies into multiple sets (recommended)
% 10 sets of frequencies for parallel computation 
% More sets can be made (number of sets is equal to the number of multiple workers available)
ind_1_eachset = floor(1: (length(freq_axis)/10) : length(freq_axis));
freq_set1 = 1 : ind_1_eachset(2)-1;
freq_set2 = ind_1_eachset(2) : ind_1_eachset(3)-1;
freq_set3 = ind_1_eachset(3) : ind_1_eachset(4)-1;
freq_set4 = ind_1_eachset(4) : ind_1_eachset(5)-1;
freq_set5 = ind_1_eachset(5) : ind_1_eachset(6)-1;
freq_set6 = ind_1_eachset(6) : ind_1_eachset(7)-1;
freq_set7 = ind_1_eachset(7) : ind_1_eachset(8)-1;
freq_set8 = ind_1_eachset(8) : ind_1_eachset(9)-1;
freq_set9 = ind_1_eachset(9) : ind_1_eachset(10)-1;
freq_set10 = ind_1_eachset(10) : length(freq_axis);

% Initialize recovered data volumes
D1_recov = zeros(nfreq,nr,ns);
D2_recov = zeros(nfreq,nr,ns);
D1_shiftrecov = zeros(nfreq,nr,ns);
D2_shiftrecov = zeros(nfreq,nr,ns);

% Initialize SNR vectors
%SNR_D1 = zeros(1,nfreq);
%SNR_D2 = zeros(1,nfreq);
%SNR_shift_D1 = zeros(1,nfreq);
%SNR_shift_D2 = zeros(1,nfreq);

%for j = freq_set1(1):freq_set1(end)  % for running multiple scripts in parallel (recommended)
for j = freq_set1(1):freq_set10(end)  % serial computation

    params.omega = 2*pi*freq_axis(j);
    opA = @(x,mode) afun(x,C,params,mode); 
    b = squeeze(D_blend_FRS(j,:,:));    
    [output,r,g,info] = spgl1(opA, b(:), 0, 0, [], opts);

    % Recover data without adding back time delays
    S1est = reshape(C'*output(1:params.curvlength),nr,ns);
    S2est = reshape(C'*output(params.curvlength+1:end),nr,ns);
    D1_recov(j,:,:) = S1est;
    D2_recov(j,:,:) = S2est;

%    SNR_D1(j) = -20*log10(norm(squeeze(D1_FRS(j,:,:)) - S1est,'fro')/norm(squeeze(D1_FRS(j,:,:)),'fro'));
%    SNR_D2(j) = -20*log10(norm(squeeze(D2_FRS(j,:,:)) - S2est,'fro')/norm(squeeze(D2_FRS(j,:,:)),'fro'));

    % Add back time delays
    S1est = exp(-1i*params.omega*params.tdelay_S1).*S1est;
    S2est = exp(-1i*params.omega*params.tdelay_S2).*S2est;
    D1_shiftrecov(j,:,:) = S1est;
    D2_shiftrecov(j,:,:) = S2est;
   
%    SNR_shift_D1(j) = -20*log10(norm(squeeze(D1_shift_FRS(j,:,:)) - S1est,'fro')/norm(squeeze(D1_shift_FRS(j,:,:)),'fro')); 
%    SNR_shift_D2(j) = -20*log10(norm(squeeze(D2_shift_FRS(j,:,:)) - S2est,'fro')/norm(squeeze(D2_shift_FRS(j,:,:)),'fro'));

end


% Save data
%rsf_write_all(fname_SNR_shift_D1, {'out=stdout'}, SNR_shift_D1)
%rsf_write_all(fname_SNR_shift_D2, {'out=stdout'}, SNR_shift_D2)
freq_st = freq_axis(freq_set1(1));  % starting frequency
rsf_write_all(fname_D1recov_FRS, {'out=stdout'}, D1_recov, [dfreq dr ds], [freq_st 0 0], {'Frequency' 'Position' 'Position'}, {'Hz' 'm' 'm'})
rsf_write_all(fname_D2recov_FRS, {'out=stdout'}, D2_recov, [dfreq dr ds], [freq_st 0 0], {'Frequency' 'Position' 'Position'}, {'Hz' 'm' 'm'})
rsf_write_all(fname_D1shiftrecov_FRS, {'out=stdout'}, D1_shiftrecov, [dfreq dr ds], [freq_st 0 0], {'Frequency' 'Position' 'Position'}, {'Hz' 'm' 'm'})
rsf_write_all(fname_D2shiftrecov_FRS, {'out=stdout'}, D2_shiftrecov, [dfreq dr ds], [freq_st 0 0], {'Frequency' 'Position' 'Position'}, {'Hz' 'm' 'm'})


% Return to the examples directory
cd(curdir)

