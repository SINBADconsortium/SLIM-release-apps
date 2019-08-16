function [opMig, opSub, dm_aj, norm_opMig] = est_opMig_norm(vel_bg, data, model_para, source, options)
% Syntax:
% [opMig,dm_aj,norm_opMig] = est_opMig_norm(vel_bg, data, model_para, source, options)
% 
% Descriptions:
% build the Born operator and roughly estimate its operator norm, can also
% output a single RTM result that is a by-product when estimating the operator
% norm. Also builds the subsampling operator.
% 
% Input list:
% vel_bg: the background velocity
% data: data you migrate (for estimating the operator norm)
% model_para: a struct that contains some modelling parameters, including grid 
%		number/distance/frequencies...
% source: a struct that contains Q/P/QmP...basically all source terms
% options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
%
% Output list:
% opMig: migration operator (Born operator)
% norm_opMig: optional. the operator norm of opMig
% dm_aj: a single RTM result
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

% acquisition parameters
nf = length(model_para.freq);
nrec = length(model_para.xrec);
nshot = length(model_para.xsrc_is);
ns_grid = length(model_para.xsrc);
nss = size(model_para.source_encoding_mat, 2);

% frequency subsampling
fidx = model_para.fidx;
nf_full = length(fidx);
opFreqSub = opRestriction(nf_full,fidx);

% source subsampling
opSourceSub = opRightMulMat([ns_grid,nshot], model_para.source_encoding_mat);

% subsampling operator
opSub = oppKron2Lo(opFreqSub, opSourceSub);

switch options.source_type
  case 1			% impulse source
    disp('Making Born operator for Greens function.')
    sourceI = opSub*source.I;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceI);
    
  case 2			% point source, wavelet estimated by EPSI
    disp('Making Born operator for primaries.')
    sourceQ = opSub*source.Q;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQ);
    
  case 3			% areal source, total wavefield
    disp('Making Born operator for multiples.')
    sourceP = opSub*source.P;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, -sourceP);
    
  case 4			% point and areal sources: wavelet by EPSI and total data
    disp('Making Born operator for total data.')
    sourceQ = opSub*source.Q;
    sourceP = opSub*source.P;
    sourceQmP = sourceQ-sourceP;
    Q = sourcegrid(ns_grid, nss, options.src_polarity, sourceQmP);
end
opMig = oppDF_old(m_bg(:), Q, model_para);

if nargout > 2
    data = opSub*data;
    dm_aj = opMig'*data;
    dm_aj = reshape(dm_aj, model_para.n);
    if nargout > 3
        data_scale = opMig*dm_aj(:);
        norm_opMig = sqrt(norm(data_scale,1)/norm(data,1));
    end
end