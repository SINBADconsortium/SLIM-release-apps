function freq = choose_freq_subset(freq_full, nb_freq)
% Syntax:
% freq = choose_freq_subset(freq_full, nb_freq)
%
% Description:
% randomly choose "nb_freq" number of frequencies from "freq_full"
%
% Input list:
% freq_full: full range of frequencies
% nb_freq: number of frequencies to choose
%
% Output list
% freq: frequencies chosen
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

% full frequency range
freq_full = vec(freq_full);
nf_full = size(freq_full,1);
% choose random subset
idx = randperm(nf_full);
idx = sort(idx(1:nb_freq));
freq = freq_full(idx); 