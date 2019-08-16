function [output,K] = timetofreq(input,ntime,nrow,ncol)

% This function performs time to frequency conversion
%
% use:
%   [output,K] = timetofreq(input,ntime,nrow,ncol)
%
% input:
% input     - input data with missing shots/receivers.
% ntime     - number of time samples in 3D data matrix
% nrow      - numbers of rows in 3D data matrix
% ncol      - numbers of columns in 3D data matrix
% output:
%  output   - frequency domain data

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


if mod(nrow,2)==0  
    output = fft(input);
    K = size(input);
else % pad the data with zeros to make source-receiver axis even
    D_new = zeros(ntime,nrow+1,ncol+1);
    D_new(:,1:nrow,1:ncol) = input;
    input = D_new;
    clear D_new
    K = size(input);
    output = fft(input);
end


% Interpolation is done only on one side of FX spectra, since other half is
% symmetric
if mod(K(1),2)==0 
    nf = floor((K(1)/2)+1);
else
    nf = ceil(K(1)/2);
end

output = output(1:nf,:,:);
K = size(output);
output = permute(reshape(output,K(1),K(2),K(3)),[3 2 1]);

end

