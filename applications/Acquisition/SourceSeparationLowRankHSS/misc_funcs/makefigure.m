function makefigure(fname, dtype, index, cax, cmap, nt, dt, nr, ns)

%------------------------------------------------------------------
% makefigure generates figures for display
%
% Use:
%   makefigure(fname, dtype, index, cax, cmap, nt, dt, nr, ns)
%
% Input:   
%     fname - name of the file with .rsf or .su extension
%     dtype - type of data to be displayed
%             (e.g. data cube in FRS or TRS domain)
%     index - time, receiver or source gather to be displayed
%       cax - color axis scale for FRS and/or TRS domain data
%      cmap - color map (default is gray)
%        nt - number of time samples
%        dt - time sampling interval
%        nr - number of receivers
%        ns - number of shots

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: June, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------


% Read data
if strfind(fname,'.rsf')
   data = rsf_read_all(fname);
elseif strfind(fname,'.su')
   data = reshape(ReadSuFast(fname,'b'), nt, nr, ns);
end

% Plot
if strcmp(dtype, 'frequency')
   figure; imagesc(squeeze(real(data(index,:,:))), [-1 1]*cax); colormap(cmap)   
   xlabel('Shot (#)'); ylabel('Receiver (#)')
elseif strcmp(dtype, 'time')
   figure; imagesc(squeeze(real(data(:,:,index))), [-1 1]*cax); colormap(cmap)      
   xlabel('Trace number'); ylabel('Time (s)')   
   ytick = 250:250:nt; 
   yticklabel = ytick*dt;
   set(gca, 'plotboxaspectratio', [1.5 2 1.5], 'Ytick', ytick, 'YTickLabel', yticklabel);
end

end  % function

