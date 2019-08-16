function [output] = syndata(input,para)

% This function will make the data for testing with missing tarces and /or noise 
% 
% use:
%   [output] = syndata(input,para)
%
% input:
% input     - fully sampled input data.
% para      - This will contain the information about data.
% output:
% output    - data with missing tarce and /or noise for testing

% Author: Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

K      = size(input);
ntime  = K(1);
nrow   = K(2);
ncol   = K(3);

if para.penalty == 1 % students't penalty function
    % Removing 50% of data
    inds           = randperm(ncol);
    ind            = inds(1:floor(ncol/para.factor)); 
    R1             = opRestriction(ncol,ind);
    R2             = opKron(R1,opDirac(nrow));
    R3             = opKron(R2,opDirac(ntime));
    input          = R3*vec(input);
    input          = reshape(input,ntime,nrow,length(ind));
    K              = size(input);
    ncolsub        = K(3);
    inds           = randperm(ncolsub);
    J              = floor(para.noise*ncolsub); % replace data by noise
    IT             = inds(1:J);
    input(:,:,IT)  = para.noiselevel*rand(ntime,nrow,length(IT)); % replace original traces with noise
    input          = R3'*vec(input);
    output         = reshape(input,ntime,nrow,ncol);

else %least-squares penalty function
    inds           = randperm(nrow);
    ind            = inds(1:floor(nrow/para.factor)); % removing 50% of data
    R1             = opRestriction(nrow,ind);
    R2             = opKron(R1,opDirac(nrow));
    R3             = opKron(R2,opDirac(ntime));
    input          = R3'*(R3*vec(input));
    output         = reshape(input,ntime,nrow,ncol);
end

end

