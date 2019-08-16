function model_data(vel, model_para, source, options)
% Syntax:
% model_data(vel, model_para, source, options)
%
% Note:
% This function is called by Migration_with_SRM. Not designed as a stand-alone
% function.
% 
% Description:
% Modelling data F(vel,Q)
%
% Input list:
% vel: the velocity model
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

% model
vel_hom = vel(1)*ones(size(vel));
m = 1./vel.^2;
m_hom = 1./vel_hom.^2;

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

% save to files
folder_name = options.results_folder;
if isempty(options.model_data_name)
    model_data_name = 'model_data.mat';
else
    model_data_name = options.model_data_name;
end

% make source
switch options.source_type
  case 1
    disp('Modelling Greens function.')
    sourceI = opSub*source.I;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceI);
  case 2
    disp('Modelling primaries.')
    sourceQ = opSub*source.Q;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQ);
  case 3
    disp('Modelling multiples.')
    sourceP = opSub*source.P;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, -sourceP);
  case 4
    disp('Modelling total data.')
    sourceQ = opSub*source.Q;
    sourceP = opSub*source.P;
    sourceQmP = sourceQ-sourceP;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQmP);
end

opFT = opKron(opDirac(nrec*nss),opFFTsym_conv_mask(options.nt, model_para.fmask));
opIFT = opFT';

% make data
D = F_old(m(:), Q, model_para);
D = reshape(gather(D), nrec, nss, nf);
D = shiftdim(D,2);
D_t = reshape(opIFT*D(:), options.nt, nrec, nss);
save([folder_name,model_data_name],'D_t');

% direct wave
if options.model_direct_wave
    D_dir = F_old(m_hom(:), Q, model_para);
    D_dir = reshape(gather(D_dir), nrec, nss, nf);
    D_dir = shiftdim(D_dir,2);
    D_dir_t = reshape(opIFT*D_dir(:), options.nt, nrec, nss);
    save([folder_name,model_data_name],'D_dir_t','-append');
end
