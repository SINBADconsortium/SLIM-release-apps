% percentage of data to contaminate
p = 0; 

% function handles to misfits used for residual and source estimation
fh1 = @(x) twonorms(x);
fh2 = @(x) twonorms(x);

% experiment label
label = 'mbase';

robust_fwi2d_marm(p,fh1,fh2,label);