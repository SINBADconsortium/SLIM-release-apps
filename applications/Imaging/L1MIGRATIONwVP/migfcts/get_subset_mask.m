function mask = get_subset_mask(quantity_full, quantity_subset, precision)
% Syntax:
% mask = get_subset_mask(quantity_full, quantity_subset)
%
% Description:
% compute mask of "quantity_subset" w.r.t "quantity_full"
%
% Input list:
% quantity_full: full range of the quantity
% quantity_subset: subset of the quantity
% precision: to what precision is the mask computed
%
% Output list:
% mask: mask, a column vector
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

if not(exist('precision','var'))
    precision = 1e-4;
end

quantity_full = round(quantity_full/precision)*precision;
quantity_subset = round(quantity_subset/precision)*precision;

% vectorize quantities
quantity_full = vec(quantity_full);
quantity_subset = vec(quantity_subset);
[~,absent_idx] = setdiff(quantity_full,quantity_subset);
mask = true(size(quantity_full));
mask(absent_idx) = false;
