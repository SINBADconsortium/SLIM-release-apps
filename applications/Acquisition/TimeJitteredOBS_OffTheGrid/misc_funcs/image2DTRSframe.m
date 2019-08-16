function image2DTRSframe(data, frame, nt, dt, caxmin, caxmax, cmap)

%----------------------------------------------------------------------------------------------------
% image2DTRSframe plots the time (T), receiver (R), and source (S) frame/slice for a given (2D) data
%
% Use:
%   image2DTRSframe(fname, frame, nt, dt, caxmin, caxmax, cmap)
%
% Input:   
%      fname - file name of a 2D data volume (with T, R, S axes) 
%      frame - a structure array including the frame/slice number to be plotted
%              (e.g., frame.t = 100, frame.r = 50, frame.s = 120)
%         nt - number of time samples
%         dt - sampling interval
%     caxmin - pseudocolor axis scaling parameter, default is minimum value of data            
%     caxmax - pseudocolor axis scaling parameter, default is maximum value of data            
%       cmap - colormap, default is gray

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2013

% Change log:
% 15 Jan 14: removed '-' (negative) sign before caxmin


% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%----------------------------------------------------------------------------------------------------

% Read the data
%data = rsf_read_all(fname);

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
time_tick_interval = 250;
time_axis = 1 : time_tick_interval : nt;
time_axis_label = (time_axis - 1)*dt;

% Time slice
figure; imagesc(real(squeeze(data(frame.t,:,:))), [caxmin caxmax]); colormap(cmap); 
set(gca,'plotboxaspectratio',[1 1 1]); colorbar;
xlabel('Shot number'); ylabel('Receiver number'); title('Time slice')

% Receiver gather
figure; imagesc(real(squeeze(data(:,frame.r,:))), [caxmin caxmax]); colormap(cmap); 
set(gca,'plotboxaspectratio',[1.5 2 1.5]); colorbar;
xlabel('Shot number'); ylabel('Time (s)'); title('Receiver gather')
set(gca,'YTick', time_axis, 'YTickLabel', time_axis_label)

% Shot gather
figure; imagesc(real(squeeze(data(:,:,frame.s))), [caxmin caxmax]); colormap(cmap); 
set(gca,'plotboxaspectratio',[1.5 2 1.5]); colorbar;
xlabel('Receiver number'); ylabel('Time (s)'); title('Shot gather')
set(gca,'YTick', time_axis, 'YTickLabel', time_axis_label)

end   % function end

