function write_ACQdata(fname, data, dt, dr, ds)

%------------------------------------------------------------------
% write_ACQdata writes output to file names defined by the user
%
% Use:
%   write_ACQdata(fname, data, dt, dr, ds)
%
% Input:   
%     fname - name of the file with .rsf extension
%      data - data to be written
%        dt - time sampling interval
%        dr - receiver sampling interval
%        ds - source sampling interval

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
%------------------------------------------------------------------

if length(size(data)) == 2
   rsf_write_all(fname, {'out=stdout'}, data, [dt dr], [0 0], {'Recording time' 'Receiver position'}, {'s' 'm'})
elseif length(size(data)) == 3
   rsf_write_all(fname, {'out=stdout'}, data, [dt dr ds], [0 0 0], {'Time' 'Receiver' 'Source'}, {'s' 'm' 'm'})
end
      
end  % function end

