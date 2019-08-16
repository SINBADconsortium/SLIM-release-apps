function [output, dim] = tfdomain(data, ntime, nfreq, nrow, ncol, mode)

%--------------------------------------------------------------------
% tfdomain performs time to frequency conversion and vice-versa.
%
% Use:
%   [output, dim] = tfdomain(data, ntime, nfreq, nrow, ncol, mode)
%
% Input:
%      data - input blended data 
%     ntime - number of time samples in 3-D data matrix
%     nfreq - number of positive frequencies in 3-D data matrix
%      nrow - numbers of rows in 3-D data matrix
%      ncol - numbers of columns in 3-D data matrix
%      mode - [1 : time to frequency], 
%             [-1 : frequency to time]
%
% Output:
%      output - frequency/time domain data
%         dim - size of output

% Author: Haneet Wason and Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: June, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%--------------------------------------------------------------------


if mode == 1  % TIME TO FREQUENCY  

   output = fft(data);
   dim = size(data);

   % Use only one side of the FX spectra, since the other half is symmetric
   if mod(dim(1),2) == 0 
      nfreq = floor((dim(1)/2)+1);
   else
      nfreq = ceil(dim(1)/2);
   end

   output = output(1:nfreq,:,:);
   dim = size(output);
   output = permute(reshape(output,dim(1),dim(2),dim(3)),[2 3 1]);

else   % FREQUENCY TO TIME
   
   if mod(ntime,2) == 0 
      output = zeros(nrow,ncol,ntime);
      output(:,:,1:nfreq) = data;
      output(:,:,1) = real(data(:,:,1));
      output(:,:,end:-1:nfreq+1) = conj(data(:,:,2:1:nfreq-1));
      output(:,:,nfreq) = real(data(:,:,nfreq));
      output = permute(output,[3 1 2]);
      output = ifft(output);
   else
      output = zeros(nrow,ncol,ntime);
      output(:,:,1:nfreq) = data;
      output(:,:,1) = real(data(:,:,1));
      output(:,:,end:-1:nfreq+1) = conj(data(:,:,2:1:nfreq));
      output = permute(output,[3 1 2]);
      output = ifft(output);
   end
   
   dim = size(output);

end % if-loop

end % function

