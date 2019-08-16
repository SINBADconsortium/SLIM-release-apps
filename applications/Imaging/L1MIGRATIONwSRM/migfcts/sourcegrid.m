function Q = sourcegrid(ns_grid, nss, src_polarity, source)
% Syntax:
% Q = sourcegrid(xs_grid, xs, src_polarity, source, source_type)
%
% Description:
% make input source for forward modeling operator F and the Born operator DF.
% It supports monopole and dipole sources. Source can be point sources or areal
% sources.
%
% Input list:
% ns_grid: number of source grid points
% nss: number of encoded sources
% src_polarity: can be 'di' or 'mo', standing for dipole and monopole sources
% source: source frequency spectrum, is a distributed vector
%
% Output list
% Q: the source.
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

% check number of inputs
if nargin < 4
    error('Not enough input arguments.')
end

% impulse response
if strcmp(src_polarity,'di')
    % dipole source
    Q_slice = zeros(2*ns_grid,ns_grid);                   
    Q_slice(1,1) = -1;
    Q_slice(2,1) = 1;
    for count = 2:ns_grid
        Q_slice(:,count) = circshift(Q_slice(:,count-1),2);
    end
else
    % monopole source
    Q_slice = eye(ns_grid);
end

% check whether source is a vector
if not(isequal(numel(source), size(source, 1)))
    error('Fatal: variable source should be a distributed vector.');
end
nf = numel(source)/(ns_grid*nss);

spmd
    source_loc = getLocalPart(source);
    nf_loc = numel(source_loc)/(ns_grid*nss);
    source_loc = reshape(source_loc, [ns_grid, nss, nf_loc]);
    % initialize Q_loc
    Q_loc = zeros(size(Q_slice,1), nss, nf_loc);
    for freq_count = 1:nf_loc
        Q_loc(:,:,freq_count) = Q_slice*source_loc(:,:,freq_count);
    end
    codistor_Q = codistributor1d(3, codistributor1d.unsetPartition, [size(Q_slice,1), nss, nf]);
    Q = codistributed.build(Q_loc, codistor_Q, 'noCommunication');
end