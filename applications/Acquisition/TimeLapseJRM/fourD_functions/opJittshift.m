function op = opJittshift(dim, jitdim, jitacq, wn)

%-------------------------------------------------------------------------------------------------------------------------------------------
% opJittshift generates the sampling operator for time-jittered OBC acquisition, for one source vessel with two airgun arrays
%
% Use:
%   opJittshift(dim, jitdim, jitacq, wn)   
%
% Input:   
%        dim - dimensions of data acquired by conventional acquisition ([nt nr ns])
%     jitdim - dimensions of data acquired by time-jittered acquisition 
%     jitacq - a structure array including the jittered acquisition parameters (output of the function: jitter_airgunarrays4D)
%         wn - frequency vector (positive and negative values)

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
gridsjitb1arr1IND = jitacq.gridsjitb1arr1IND;
gridsjitb1arr2IND = jitacq.gridsjitb1arr2IND;
tfirejitb1arr1gridIND = jitacq.tfirejitb1arr1gridIND;
tfirejitb1arr2gridIND = jitacq.tfirejitb1arr2gridIND;
tshiftb1arr1 = jitacq.tshiftb1arr1;
tshiftb1arr2 = jitacq.tshiftb1arr2;

% Sampling operator
fh = @(x, mode) opJittshift_intrnl(m, n, x, mode, wn, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2);
op = opFunction(m, n, fh);


%=========================================================================%


function y = opJittshift_intrnl(m, n, x, mode, wn, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2)

Ft = opDFT(dim(1));

if (mode == 0)
    y = {m, n, [0,0,0,0], {'opJittshift'}};   
elseif (mode == 1)
    x = reshape(x,dim);
    y1 = zeros(jitdim);
    y2 = y1;
    
    for j = 1 : length(gridsjitb1arr1IND)
        shot_arr1 = x(:, :, gridsjitb1arr1IND(j));
        shot_arr1 = Ft*shot_arr1;
        for k = 1:size(shot_arr1,2)
	    shot_arr1(:,k) = shot_arr1(:,k).*exp(-1i*wn*tshiftb1arr1(j));
        end
        shot_arr1 = Ft'*shot_arr1;

        shot_arr2 = x(:, :, gridsjitb1arr2IND(j));
        shot_arr2 = Ft*shot_arr2;
        for k = 1:size(shot_arr2,2)
            shot_arr2(:,k) = shot_arr2(:,k).*exp(-1i*wn*tshiftb1arr2(j));
        end
        shot_arr2 = Ft'*shot_arr2;

        y1(tfirejitb1arr1gridIND(j) : tfirejitb1arr1gridIND(j) + dim(1) - 1, :) = y1(tfirejitb1arr1gridIND(j) : tfirejitb1arr1gridIND(j) + dim(1) - 1, :) + shot_arr1;
        y2(tfirejitb1arr2gridIND(j) : tfirejitb1arr2gridIND(j) + dim(1) - 1, :) = y2(tfirejitb1arr2gridIND(j) : tfirejitb1arr2gridIND(j) + dim(1) - 1, :) + shot_arr2;
    end
     
    y = y1 + y2;
    y = y(:);   
else
    x1 = reshape(x, jitdim);
    x2 = x1;
    y1 = zeros(dim);
    y2 = y1;
    
    for j = 1 : length(gridsjitb1arr1IND)     
      y1(:,:,gridsjitb1arr1IND(j)) = x1(tfirejitb1arr1gridIND(j) : tfirejitb1arr1gridIND(j) + dim(1) - 1, :);
      shot_arr1 = Ft*y1(:,:,gridsjitb1arr1IND(j));
      for k = 1:size(shot_arr1,2)
          shot_arr1(:,k) = shot_arr1(:,k).*exp(-1i*wn*-tshiftb1arr1(j));
      end   
      y1(:,:,gridsjitb1arr1IND(j)) = Ft'*shot_arr1; 
      
      y2(:,:,gridsjitb1arr2IND(j)) = x2(tfirejitb1arr2gridIND(j) : tfirejitb1arr2gridIND(j) + dim(1) - 1, :);
      shot_arr2 = Ft*y2(:,:,gridsjitb1arr2IND(j));
      for k = 1:size(shot_arr2,2)
          shot_arr2(:,k) = shot_arr2(:,k).*exp(-1i*wn*-tshiftb1arr2(j));
      end   
      y2(:,:,gridsjitb1arr2IND(j)) = Ft'*shot_arr2; 
    end 

    y = y1 + y2;
    y = y(:);  
end

end  % internal function end

end  % external function end

