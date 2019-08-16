function save_data(dtype, fname, data, samp_int)

%------------------------------------------------------------------
% save_data saves the output to file names defined by the user
%
% Use:
%   save_data(dtype, fname, data, samp_int)
%
% Input:   
%        dtype - type of data to be saved 
%                (e.g. data cube in FRS or TRS domain)
%        fname - name of the file with .rsf extension
%         data - data to be written
%     samp_int - frequency or time sampling interval

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------

if strcmp(dtype, 'data_FRS')
   rsf_write_all(fname, {'out=stdout'}, data, [samp_int 1 1], [0 0 0], {'Frequency' 'Trace' 'Trace'}, {'Hz' '#' '#'})
elseif strcmp(dtype, 'data_TRS')
   rsf_write_all(fname, {'out=stdout'}, data, [samp_int 1 1], [0 0 0], {'Time' 'Trace' 'Trace'}, {'s' '#' '#'})
end
      
end  % function end

