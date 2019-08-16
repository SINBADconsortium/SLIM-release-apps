function migvp_testop
% Syntax:
% test_migsrm_fcts
%
% Descriptions:
% test operators and functions used for the algorithm
%
% No input and output
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Apr/25/2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% test get_subset_mask
dt = 0.004;
nt = 256;
df = 1/(nt*dt);
fu = 1/(2*dt);
freq_full = 0:df:fu;
nf = length(freq_full);
fmask_true = false(nf,1);
idx = randperm(nf);
fmask_true(idx(1:10)) = true;
freq = freq_full(fmask_true);
fmask_calc = get_subset_mask(freq_full,freq);
if isequal(fmask_true,fmask_calc)
	disp('Test passed for masking function.')
else
	disp('Test failed for masking function.')
end
clear

% test opRealRestriction
opRealRes = opRealRestriction(32);
flag_dottest = not(dottest(opRealRes,10));

test_var = rand(32,1)+1i*rand(32,1);
out_var = opRealRes*test_var;
flag_unit = isreal(out_var);
if flag_dottest && flag_unit
	disp('Test passed for opRealRestriction.')
else
	disp('Test failed for opRealRestriction.')
end
clear

% test opRightMulMat
W = rand(32,5);
opEncode = opRightMulMat([32,32],W);
flag_dottest = not(dottest(opEncode,10));

test_var = rand(32,32);
out_var = opEncode*test_var(:);
out_var = reshape(out_var,32,5);
flag_unit = isequal(out_var,test_var*W);
if flag_dottest && flag_unit
	disp('Test passed for opRightMulMat.')
else
	disp('Test failed for opRightMulMat.')
end
clear

% opTriMute is by Tim Lin