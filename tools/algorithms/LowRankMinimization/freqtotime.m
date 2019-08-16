function [ Dt,K ] = freqtotime(data,para,ntime)

% This function performs frequency to time domain conversion
%
% use:
%   [ Drecover,K ] = freqtotime(data,para,ntime)
%
% input:
% data      - input data in frequency domain
% para      - This will contain the information about data.
% ntime     - number of time samples in 3D data matrix
% output:
%  Dt - time domain data

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

if mod(ntime,2)==0 
    Dt = zeros(para.nrow,para.ncol,ntime);
    Dt(:,:,1:para.nf)= data;
    Dt(:,:,1)=real(data(:,:,1));
    Dt(:,:,end:-1:para.nf+1)=conj(data(:,:,2:1:para.nf-1));
    Dt(:,:,para.nf)=real(data(:,:,para.nf));
    Dt = permute(Dt,[3 2 1]);
    Dt = ifft(Dt);
    if para.odd==1 
    Dt = Dt(:,1:para.nrow-1,1:para.ncol-1);
    end
else
    Dt = zeros(para.nrow,para.ncol,ntime);
    Dt(:,:,1:para.nf)= data;
    Dt(:,:,1)=real(data(:,:,1));
    Dt(:,:,end:-1:para.nf+1)=conj(data(:,:,2:1:para.nf));
    Dt = permute(Dt,[3 2 1]);
    Dt = ifft(Dt);
    if para.odd==1  
    Dt = Dt(:,1:para.nrow-1,1:para.ncol-1);
    end 
end

K = size(Dt);

end

