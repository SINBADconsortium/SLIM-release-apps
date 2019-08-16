function test_model(vel_bg, vel_true, model_para, source, options)
% Syntax:
% test_model(model_bg, model_true, model_para, source, options)
%
% Note:
% This function is called by Migration_with_SRM. Not designed as a stand-alone
% function.
% 
% Description:
% Test the following.
% 1. whether the background velocity/density models are smooth enough not to give rise to any reflection
% 2. whether the absorbing boundary is wide enough for wave at boundary to die out
% 3. difference between the linearized data and the nonlinear residual [F(m)-F(m0)]
%
% Input list:
% vel_bg: the background velocity
% vel_true: true velocity
% model_para: a struct that contains some modelling parameters, including grid 
%		number/distance/frequencies...
% source: a struct that contains source.I, source.Q and source.P
% options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
%
% Output List:
% results will be saved to file.
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

% models
m_bg = 1./vel_bg.^2;
m_true = 1./vel_true.^2;
m_pert = m_true - m_bg;

% acquisition parameters
nf = length(model_para.freq);
nrec = length(model_para.xrec);
nshot = length(model_para.xsrc_is);
ns_grid = length(model_para.xsrc);
nss = size(model_para.source_encoding_mat, 2);

% sampling over frequencies
fidx = model_para.fidx;
nf_full = length(fidx);
opFreqSub = opRestriction(nf_full,fidx);

% sampling over sources
opSourceSub = opRightMulMat([ns_grid,nshot], model_para.source_encoding_mat);

% subsampling operator
opSub = oppKron2Lo(opFreqSub, opSourceSub);

% make source
switch options.source_type
  case 1
    disp('Test using Greens function.')
    sourceI = opSub*source.I;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceI);
  case 2
    disp('Test using primaries.')
    sourceQ = opSub*source.Q;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQ);
  case 3
    disp('Test using multiples.')
    sourceP = opSub*source.P;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, -sourceP);
  case 4
    disp('Test using total data.')
    sourceQ = opSub*source.Q;
    sourceP = opSub*source.P;
    sourceQmP = sourceQ-sourceP;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQmP);
end

opFT = opKron(opDirac(nrec*nss),opFFTsym_conv_mask(options.nt, model_para.fmask));
opIFT = opFT';

% first test
[D_bg, J] = F_old(m_bg(:), Q, model_para);
D_bg = reshape(gather(D_bg), nrec, nss, nf);
D_bg = shiftdim(D_bg,2);
D_bg_t = reshape(opIFT*D_bg(:), options.nt, nrec, nss);

% second test
D_true = F_old(m_true(:), Q, model_para);
D_true = reshape(gather(D_true), nrec, nss, nf);
D_true = shiftdim(D_true,2);
D_true_t = reshape(opIFT*D_true(:), options.nt, nrec, nss);

% third test
D_linear = J*m_pert(:);
D_linear = reshape(gather(D_linear), nrec, nss, nf);
D_linear = shiftdim(D_linear,2);
D_linear_t = reshape(opIFT*D_linear(:), options.nt, nrec, nss);

%% save test files
folder_name = options.results_folder;
if isempty(options.test_file_name)
    test_file_name = 'test_model_data.mat';
else
    test_file_name = options.test_file_name;
end
save([folder_name,test_file_name],'D_bg_t','D_true_t','D_linear_t');