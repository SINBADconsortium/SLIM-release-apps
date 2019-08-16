function make_linear_data(vel_bg, m_pert, model_para, source, options)
% Syntax:
% make_linear_data(vel_bg, m_pert, model_para, source, options)
%
% Note:
% This function is called by Migration_with_SRM. Not designed as a stand-alone
% function.
%
% Description:
% Make linear data: Green's function/primaries/multiples/total data depending on
% your input data type.
% 1.) data_type = 1: make Green's function
% 2.) data_type = 2: make primary
% 3.) data_type = 3: make multiples
% 4.) data_type = 4: make total data
%
% Input list:
% vel_bg: the background velocity
% m_pert: model perturbation
% model_para: a struct that contains some modelling parameters, including grid 
%		number/distance/frequencies...
% source: a struct that contains Q/P/QmP...basically all source terms
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

% i/o
folder_name = options.results_folder;
if isempty(options.linear_data_name)
    linear_data_name = 'linear_data.mat';
else
    linear_data_name = options.linear_data_name;
end

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

opFT = opKron(opDirac(nrec*nss),opFFTsym_conv_mask(options.nt, model_para.fmask));
opIFT = opFT';

switch options.data_type
  case 1			% impulse source
    disp('Making linearized Greens function.')
    sourceI = opSub*source.I;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceI);
    opMig = oppDF_old(m_bg(:), Q, model_para);
    GF_linear_F = gather(opMig*m_pert(:));
    GF_linear_F = reshape(GF_linear_F, nrec, nss, nf);
    GF_linear_F = shiftdim(GF_linear_F,2);
    GF_linear = reshape(opIFT*GF_linear_F(:), options.nt, nrec, nss);
    
    if ~exist([folder_name,linear_data_name],'file')
        save([folder_name,linear_data_name],'GF_linear','GF_linear_F')
    else
        save([folder_name,linear_data_name],'GF_linear','GF_linear_F','-append')
    end
    
  case 2			% point source, wavelet estimated by EPSI
    disp('Making linearized primaries.')
    sourceQ = opSub*source.Q;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQ);
    opMig = oppDF_old(m_bg(:), Q, model_para);
    primary_linear_F = gather(opMig*m_pert(:));
    primary_linear_F = reshape(primary_linear_F, nrec, nss, nf);
    primary_linear_F = shiftdim(primary_linear_F,2);
    primary_linear = reshape(opIFT*primary_linear_F(:), options.nt, nrec, nss);
    if ~exist([folder_name,linear_data_name],'file')
        save([folder_name,linear_data_name],'primary_linear','primary_linear_F')
    else
        save([folder_name,linear_data_name],'primary_linear','primary_linear_F','-append')
    end
    
  case 3			% areal source, total wavefield
    disp('Making linearized multiples.')
    sourceP = opSub*source.P;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, -sourceP);
    opMig = oppDF_old(m_bg(:), Q, model_para);
    multiple_linear_F = gather(opMig*m_pert(:));
    multiple_linear_F = reshape(multiple_linear_F, nrec, nss, nf);
    multiple_linear_F = shiftdim(multiple_linear_F,2);
    multiple_linear = reshape(opIFT*multiple_linear_F(:), options.nt, nrec, nss);
    if ~exist([folder_name,linear_data_name],'file')
        save([folder_name,linear_data_name],'multiple_linear','multiple_linear_F')
    else
        save([folder_name,linear_data_name],'multiple_linear','multiple_linear_F','-append')
    end
    
  case 4			% point and areal sources: wavelet by EPSI and total data
    disp('Making linearized total up-going data.')
    sourceQ = opSub*source.Q;
    sourceP = opSub*source.P;
    sourceQmP = sourceQ-sourceP;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQmP);
    opMig = oppDF_old(m_bg(:), Q, model_para);
    totaldata_linear_F = gather(opMig*m_pert(:));
    totaldata_linear_F = reshape(totaldata_linear_F, nrec, nss, nf);
    totaldata_linear_F = shiftdim(totaldata_linear_F,2);
    totaldata_linear = reshape(opIFT*totaldata_linear_F(:), options.nt, nrec, nss);
    if ~exist([folder_name,linear_data_name],'file')
        save([folder_name,linear_data_name],'totaldata_linear','totaldata_linear_F')
    else
        save([folder_name,linear_data_name],'totaldata_linear','totaldata_linear_F','-append')
    end
end