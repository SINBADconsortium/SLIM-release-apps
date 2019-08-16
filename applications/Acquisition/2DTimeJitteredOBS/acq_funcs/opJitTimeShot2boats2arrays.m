function op = opJitTimeShot2boats2arrays(dim, jitdim, jitacq)

%--------------------------------------------------------------------------------------------------------------------------------------------------
% opJitTimeShot2boats2arrays generates the sampling operator for time-jittered OBC acquisition, for two source vessels with two airgun arrays each
%
% Use:
%   opJitTimeShot2boats2arrays(dim, jitdim, jitacq)   
%
% Input:   
%        dim - dimensions of data acquired by conventional acquisition ([nt nr ns])
%     jitdim - dimensions of data acquired by time-jittered acquisition 
%     jitacq - a structure array including the jittered acquisition parameters (output of the function: jitter_airgunarrays)      

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2013

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%--------------------------------------------------------------------------------------------------------------------------------------------------

% Size of the sampling operator
m = prod(jitdim);
n = prod(dim);

% Jittered acquisition parameters
gridsjitb1arr1IND = jitacq.gridsjitb1arr1IND;
gridsjitb1arr2IND = jitacq.gridsjitb1arr2IND;
gridsjitb2arr1IND = jitacq.gridsjitb2arr1IND;
gridsjitb2arr2IND = jitacq.gridsjitb2arr2IND;
tfirejitb1arr1gridIND = jitacq.tfirejitb1arr1gridIND;
tfirejitb1arr2gridIND = jitacq.tfirejitb1arr2gridIND;
tfirejitb2arr1gridIND = jitacq.tfirejitb2arr1gridIND;
tfirejitb2arr2gridIND = jitacq.tfirejitb2arr2gridIND;

% Sampling operator
fh = @(x, mode) opJitTimeShot2boats2arrays_intrnl(m, n, x, mode, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, gridsjitb2arr1IND, gridsjitb2arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tfirejitb2arr1gridIND, tfirejitb2arr2gridIND);
op = opFunction(m, n, fh);


%=========================================================================%


function y = opJitTimeShot2boats2arrays_intrnl(m, n, x, mode, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, gridsjitb2arr1IND, gridsjitb2arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tfirejitb2arr1gridIND, tfirejitb2arr2gridIND)


if (mode == 0)
    y = {m, n, [0,0,0,0], {'opJitTimeShot2boats2arrays'}};   
elseif (mode == 1)
    x = reshape(x,dim);
    y1 = zeros(jitdim);
    y2 = y1;
    y3 = y1;
    y4 = y1;
    
    for j = 1 : length(gridsjitb1arr1IND)
        shot_b1arr1 = x(:, :, gridsjitb1arr1IND(j));
        shot_b1arr2 = x(:, :, gridsjitb1arr2IND(j));
        shot_b2arr1 = x(:, :, gridsjitb2arr1IND(j));
        shot_b2arr2 = x(:, :, gridsjitb2arr2IND(j));
        y1(tfirejitb1arr1gridIND(j) : tfirejitb1arr1gridIND(j) + dim(1) - 1, :) = y1(tfirejitb1arr1gridIND(j) : tfirejitb1arr1gridIND(j) + dim(1) - 1, :) + shot_b1arr1;
        y2(tfirejitb1arr2gridIND(j) : tfirejitb1arr2gridIND(j) + dim(1) - 1, :) = y2(tfirejitb1arr2gridIND(j) : tfirejitb1arr2gridIND(j) + dim(1) - 1, :) + shot_b1arr2;
        y3(tfirejitb2arr1gridIND(j) : tfirejitb2arr1gridIND(j) + dim(1) - 1, :) = y3(tfirejitb2arr1gridIND(j) : tfirejitb2arr1gridIND(j) + dim(1) - 1, :) + shot_b2arr1;
        y4(tfirejitb2arr2gridIND(j) : tfirejitb2arr2gridIND(j) + dim(1) - 1, :) = y4(tfirejitb2arr2gridIND(j) : tfirejitb2arr2gridIND(j) + dim(1) - 1, :) + shot_b2arr2;
    end
     
    y = y1 + y2 + y3 + y4;
    y = y(:);   
else
    x1 = reshape(x, jitdim);
    x2 = x1;
    x3 = x1;
    x4 = x1;
    y1 = zeros(dim);
    y2 = y1;
    y3 = y1;
    y4 = y1;    

    for j = 1 : length(gridsjitb1arr1IND)
      y1(:,:,gridsjitb1arr1IND(j)) = x1(tfirejitb1arr1gridIND(j) : tfirejitb1arr1gridIND(j) + dim(1) - 1, :);
      y2(:,:,gridsjitb1arr2IND(j)) = x2(tfirejitb1arr2gridIND(j) : tfirejitb1arr2gridIND(j) + dim(1) - 1, :);
      y3(:,:,gridsjitb2arr1IND(j)) = x3(tfirejitb2arr1gridIND(j) : tfirejitb2arr1gridIND(j) + dim(1) - 1, :);
      y4(:,:,gridsjitb2arr2IND(j)) = x4(tfirejitb2arr2gridIND(j) : tfirejitb2arr2gridIND(j) + dim(1) - 1, :);
    end 

    y = y1 + y2 + y3 + y4;
    y = y(:);  
end

end  % internal function end

end  % external function end


