function encoding_mat = make_encoding_mat(nshot, nss, type)
% Syntax:
% encoding_mat = make_encoding_mat(nshot, nss, type)
%
% Description:
% make source encoding matrix
%
% Input list:
% nshot: number of shots
% nss: number of encoded sources
%
% Output list
% encoding_mat: encoding matrix
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: May/07/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if strcmp(type,'gaussian')
    % use random gaussian distributed weighting
    encoding_mat = randn(nshot, nss);
elseif strcmp(type,'dirac')
    % choose random subset from all sequential sources
    idx = randperm(nshot);
    idx = sort(idx(1:nss));
    encoding_mat = eye(nshot);
    encoding_mat = encoding_mat(:,idx);
end