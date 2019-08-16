function Dmap = fdata2AcqGrid(inputData, mode, ncorig, startIND, nr_subset, nc_subset)

%------------------------------------------------------------------------------
% fdata2AcqGrid reads .mat or .rsf data in (offset, source) or (frequency, offset, source) coordinates and embeds it in (receiver, source) array, i.e., it maps data from dynamic (moving receivers) geometry to static (fixed receivers) geometry.
% 
% Use:
%   fdata2AcqGrid(inputData, mode, ncorig, startIND, nr_subset, nc_subset)  
%
% Input: 
%     inputData - input data in .mat or .rsf format
%          mode - 1: dynamic to static geometry, 2: static to dynamic geometry
%        ncorig - number of columns in mapped data (i.e., static geometry) 
%      startIND - column index corresponding to the first full-offset location
%     nr_subset - number of rows corresponding to data with full offset 
%                 (same as the number of rows of input data)          
%     nc_subset - number of columns corresponding to data with full offset 
%                 (same as the number of columns of input data)          
%
% Output:
%      Dmap - data mapped on to static or dynamic grid 

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
% Contact: Haneet Wason (hwason@eos.ubc.ca)
%         
% Date: April, 2012

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------------------

if ischar(inputData)
   if strfind(inputData,'.rsf')
      data = rsf_read_all(inputData);
      data = squeeze(data);
   elseif strfind(inputData,'.mat')
      load(inputData);
      data = Uobs;
      clear Uobs
   end
else
   data = inputData; 
end

if mode == 1 % convert to static geometry
    [n1,n2] = size(data);
    Dmap = zeros(n1+n2-1);
    k = 0;
    for i = 1:n1
        for j = 1:n2
            Dmap(j+n1-1,j+n1-1+k) = data(i,j);
        end
        k = k-1;
    end
    nb = size(Dmap,1);
    for i = 1:nb
        for j = 1:nb
            Dmap(i,j) = Dmap(j,i);
        end
    end    
elseif mode == 2  % convert back to dynamic geometry
    Dmap = zeros(nr_subset, nc_subset);
    q = 0;
    k = 1;
    
    for j = startIND:ncorig
        Dmap(:,k) = data(j : -1 : 1+q, j);
        k = k+1;
        q = q+1;
    end
end  

end  % function

