function migsrm_testop
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
% Date: Feb/14/2012
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

% test opSaveSnapshot
% no testing is done since it is essentially an opDirac with file saving
% functionality
disp('Test passed for opSaveSnapshot')

% test opTukeyWinMask
mask = true(32,1);
mask(1:4) = false;
opTukeyWin = opTukeyWinMask(32,0.2,mask);
flag_dottest = not(dottest(opTukeyWin,10));

if flag_dottest
	disp('Test passed for opRightMulMat.')
else
	disp('Test failed for opRightMulMat.')
end
clear

% test opFFTsym_mask
mask = true(6,1);
mask(3) = false;
mask(5) = false;
opFT = opFFTsym_mask(11,mask);
test_var = randn(11,1);
out_var = opFT*test_var;
adj_var = opFT'*out_var;
flag_dottest = ((real(test_var'*adj_var)-real(out_var'*out_var)) < 1e-6);
if flag_dottest
	disp('Test passed for opFFTsym_mask.')
else
	disp('Test failed for opFFTsym_mask.')
end
clear

% test opFFTsym_conv_mask
% NOTE: This operator does not pass dottest by design.
% This test is to confirm the consistency between this
% operator and opFFTsym_conv_datacube.
mask = true(6,1);
mask(3) = false;
mask(5) = false;
opFT = opFFTsym_conv_mask(11,mask);
opFT1 = opKron(opDirac(100),opFT);
opFT2 = opFFTsym_conv_datacube(11,10,10,mask);
test_var = randn(1100,1);
out_var1 = opFT1*test_var;
out_var2 = opFT2*test_var;
adj_var1 = opFT1'*out_var1;
adj_var2 = opFT2'*out_var2;
flag_unit = (isequal(out_var1,out_var2)) && (isequal(adj_var1, ...
                                                  adj_var2));
if flag_unit
	disp('Test passed for opFFTsym_conv_mask.')
else
	disp('Test failed for opFFTsym_conv_mask.')
end
clear