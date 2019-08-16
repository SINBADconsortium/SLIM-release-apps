% Convert recovered data from frequency domain to time domain

% Set a label for the experiment
label = 'SourceSepL1';

% Set paths
setpath;
cd(expdir);

% Number of samples
nt = 1250;
nr = 231;
ns = 231;

% Sampling intervals
dt = 0.004;
dfreq = 1/(nt*dt);
dr = 20;
ds = 20;

% Use only one side of the FX spectra, since the other half is symmetric
if mod(nt,2) == 0 
   nfreq = floor((nt/2)+1);
else
   nfreq = ceil(nt/2);
end

% Load FRS results
% FRS: data in frequency-receiver-source coordinates
fname_D1recov_FRS = [label '_src1recov_FRS.rsf'];
fname_D2recov_FRS = [label '_src2recov_FRS.rsf'];
fname_D1shiftrecov_FRS = [label '_shifted_src1recov_FRS.rsf'];
fname_D2shiftrecov_FRS = [label '_shifted_src2recov_FRS.rsf'];
D1recov_FRS = rsf_read_all(fname_D1recov_FRS);
D2recov_FRS = rsf_read_all(fname_D2recov_FRS);
D1shiftrecov_FRS = rsf_read_all(fname_D1shiftrecov_FRS);
D2shiftrecov_FRS = rsf_read_all(fname_D2shiftrecov_FRS);

% Frequency to time conversion
D1recov_TRS = tfdomain(permute(D1recov_FRS, [2 3 1]), nt, nfreq, nr, nc, -1);
D2recov_TRS = tfdomain(permute(D2recov_FRS, [2 3 1]), nt, nfreq, nr, nc, -1);
D1shiftrecov_TRS = tfdomain(permute(D1shiftrecov_FRS, [2 3 1]), nt, nfreq, nr, nc, -1);
D2shiftrecov_TRS = tfdomain(permute(D2shiftrecov_FRS, [2 3 1]), nt, nfreq, nr, nc, -1);

% Save data
% TRS: data in time-receiver-source coordinates
fname_D1recov_TRS = [label '_src1recov_TRS.rsf'];
fname_D2recov_TRS = [label '_src2recov_TRS.rsf'];
fname_D1shiftrecov_TRS = [label '_shifted_src1recov_FRS.rsf'];
fname_D2shiftrecov_TRS = [label '_shifted_src2recov_FRS.rsf'];
rsf_write_all(fname_D1recov_TRS, {'out=stdout'}, D1recov_TRS, [dt dr ds], [0 0 0], {'Time' 'Position' 'Position'}, {'s' 'm' 'm'})
rsf_write_all(fname_D2recov_TRS, {'out=stdout'}, D2recov_TRS, [dt dr ds], [0 0 0], {'Time' 'Position' 'Position'}, {'s' 'm' 'm'})
rsf_write_all(fname_D1shiftrecov_TRS, {'out=stdout'}, D1shiftrecov_TRS, [dt dr ds], [0 0 0], {'Time' 'Position' 'Position'}, {'s' 'm' 'm'})
rsf_write_all(fname_D2shiftrecov_TRS, {'out=stdout'}, D2shiftrecov_TRS, [dt dr ds], [0 0 0], {'Time' 'Position' 'Position'}, {'s' 'm' 'm'})

