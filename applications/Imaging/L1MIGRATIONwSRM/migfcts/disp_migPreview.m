function disp_migPreview(preview_file)
% Syntax
% disp_migPreview(preview_file)
%
% Description
% show results
%
% Input list:
% preview_file: full path of preview file
%
% Output list: None. Display plots.
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

load(preview_file)
plot_utils.show_model(dm, options)
title('Inverted model perturbation')