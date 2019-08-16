function imageJitteredData(fname, t1, t2, dt, nr, dr, caxmin, caxmax, cmap)

%--------------------------------------------------------------------------------------------------------------------
% imageJitteredData plots a section of the blended data, t1 - t2 seconds from the full recording time
%
% Use: 
%   imageJitteredData(fname, t1, t2, dt, nr, dr, caxmin, caxmax, cmap)
%
% Input:
%      fname - file name for jittered data volume (a 2D matrix of size: <recording time units, number of receivers>)
%         t1 - start time for the section to be plotted (in seconds)
%         t2 - end time for the section to be plotted (in seconds)
%         dt - sampling interval
%         nr - number of receivers
%         dr - receiver sampling interval
%     caxmin - pseudocolor axis scaling parameter, default is minimum value of data            
%     caxmax - pseudocolor axis scaling parameter, default is maximum value of data 
%       cmap - colormap, default is gray

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
%--------------------------------------------------------------------------------------------------------------------

% Read the data
data = rsf_read_all(fname);

if not(exist('caxmin', 'var')) 
  caxmin = min(data(:));
end 

if not(exist('caxmax', 'var')) 
  caxmax = max(data(:));
end

if not(exist('cmap', 'var')) 
   cmap = 'gray';
end

% Labeling the time axis
time_tick_interval = 5.0;
time_tick_persamp  = floor(time_tick_interval/dt);
time_axis          = 1 : time_tick_persamp : (t2-t1)/dt + 1;
time_axis_label    = (t1 : time_tick_interval : t2);

% Labeling the receiver axis
rec_tick_interval = 40;
rec_axis          = 1 : rec_tick_interval : nr;
rec_axis_label    = (rec_axis - 1)*dr;

% View the data
figure
imagesc(data(t1/dt:t2/dt, :), [-caxmin caxmax]); colormap(cmap); colorbar;
xlabel('Receiver position (m)'); ylabel('Recording time (s)'); title('Jittered data')
set(gca, 'plotboxaspectratio', [1.5 2 1.5]);
set(gca, 'YTick', time_axis, 'YTickLabel', time_axis_label, 'XTick', rec_axis, 'XTickLabel', rec_axis_label);

end  % function end

