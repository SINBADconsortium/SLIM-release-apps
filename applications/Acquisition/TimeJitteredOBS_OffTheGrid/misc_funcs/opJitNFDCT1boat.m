function op = opJitNFDCT1boat(dim, jitdim, jitacq, idx, wn, nr)

% Designed to work on data generated via NFFT

%-------------------------------------------------------------------------------------------------------------------------------------------
% opJitNFDCT1boat generates the sampling operator for time-jittered OBC acquisition, for one source vessel with two airgun arrays
%
% Use:
%   opJitNFDCT1boat(dim, jitdim, jitacq, idx, wn, nr)   
%
% Input:   
%        dim - dimensions of data acquired by conventional acquisition ([nt nr ns])
%     jitdim - dimensions of data acquired by time-jittered acquisition 
%     jitacq - a structure array including the jittered acquisition parameters (output of the function: jitter_airgunarrays)      
%        idx -
%         wn - frequency samples vector
%         nr - number of receivers

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: June, 2013

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%-------------------------------------------------------------------------------------------------------------------------------------------


% Size of the sampling operator
m = prod(jitdim);
n = prod(dim);

% Jittered acquisition parameters
tfirejitb1arr1gridIND = jitacq.tfirejitb1arr1gridIND;
tfirejitb1arr2gridIND = jitacq.tfirejitb1arr2gridIND;
tshiftb1arr1 = jitacq.tshiftb1arr1;
tshiftb1arr2 = jitacq.tshiftb1arr2;

% Sampling operator
fh = @(x, mode) opJitNFDCT1boat_intrnl(m, n, x, mode, dim, jitdim, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2, idx, wn, nr);
op = opFunction(m, n, fh);


%=========================================================================%


function y = opJitNFDCT1boat_intrnl(m, n, x, mode, dim, jitdim, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2, idx, wn, nr)

tfirejitgridIND = ([tfirejitb1arr1gridIND tfirejitb1arr2gridIND]);
tshift = ([tshiftb1arr1 tshiftb1arr2]);

%Ft = opFFT1C(dim(1));  % Ft': analysis, Ft: synthesis
Ft = opKron(opDirac(nr),opDFT(dim(1)));

if (mode == 0)

   y = {m, n, [0,0,0,0], {'opJitNFDCT1boat'}};   

elseif (mode == 1)

   x = reshape(x,dim);
   y = zeros(jitdim);

   for j = 1 : length(tfirejitgridIND)
       shot = x(:, :, j);    
       Ftshot = reshape(Ft*shot(:), size(shot));
       for k = 1:size(Ftshot,2)
           Ftshot(:,k) = Ftshot(:,k).*exp(-1i*wn*tshift(idx(j)));
       end
       x(:,:,j) = reshape(Ft'*Ftshot(:),size(shot));
       y(tfirejitgridIND(idx(j)) : tfirejitgridIND(idx(j)) + dim(1) - 1, :) = y(tfirejitgridIND(idx(j)) : tfirejitgridIND(idx(j)) + dim(1) - 1, :) + x(:,:,j);
   end

   y = y(:);   

else

   x = reshape(x, jitdim);
   y = zeros(dim);
    
   for j = 1 : length(tfirejitgridIND)
       y(:,:,j) = x(tfirejitgridIND(idx(j)) : tfirejitgridIND(idx(j)) + dim(1) - 1, :);
       shot = y(:,:,j);
       Ftshot = reshape(Ft*shot(:),size(shot));    
       for k = 1:size(Ftshot,2)
   	   Ftshot(:,k) = Ftshot(:,k).*exp(-1i*wn*-tshift(idx(j)));
       end
       y(:,:,j) = reshape(Ft'*Ftshot(:),size(shot));       
   end 

   y = y(:);  

end

end  % internal function end

end  % external function end

