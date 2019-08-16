function write_EPSIwavelet(filename,wavelet,nt,dt,wavelet_window_start,wavelet_window_end,wavelet_invObliq)
% Writes a single trace

% Check filetype
if strcmp(filename(end-3:end),'.mat')
    % matlabfile
    filetype = 'matlab';
    filename_invObliq = [filename(1:end-4) '.corrected' filename(end-3:end)];
    disp('Writing source wavelet to a matlab datafile...')
elseif strcmp(filename(end-3:end),'.bin')
    % pure binary type
    filetype = 'binary';
    filename_invObliq = [filename(1:end-4) '.corrected' filename(end-3:end)];
	disp('Writing source wavelet to a binary datafile...')
elseif strcmp(filename(end-2:end),'.su')
    filetype = 'su';
    filename_invObliq = [filename(1:end-3) '.corrected' filename(end-2:end)];
    disp('Writing source wavelet to an su datafile...')
else
    error(['Unrecognized filetype for output file: ' filename])
end

% Preprocess wavelet, window it out
wavelet = wavelet(:);
wavelet(wavelet_window_end+1:end) = [];
wavelet(1:wavelet_window_start-1) = [];

if exist('wavelet_invObliq','var')
    wavelet_invObliq = wavelet_invObliq(:);
    wavelet_invObliq(wavelet_window_end+1:end) = [];
    wavelet_invObliq(1:wavelet_window_start-1) = [];
end

% Determine start time
lagtime = (wavelet_window_start-1-nt)*dt;
endtime = (wavelet_window_end-nt-1)*dt;
time_axis = [lagtime:dt:endtime];


switch filetype
    
    case 'matlab'
        save(filename,'wavelet','dt','time_axis')
        if exist('wavelet_invObliq','var')
            wavelet = wavelet_invObliq;
            save(filename_invObliq,'wavelet','dt','time_axis')
        end
    case 'binary'
        io_utils.WriteBin(filename,wavelet(:))
        if exist('wavelet_invObliq','var')
            io_utils.WriteBin(filename_invObliq,wavelet_invObliq(:))
        end
        disp(['    lagtime = ' num2str(lagtime)])
        disp(['    dt = ' num2str(dt)])
        
    case 'su'
        io_utils.WriteSu(filename,wavelet(:),'dt',dt,'DelayRecordingTime',lagtime*1e3) % lagtime is in milliseconds
        if exist('wavelet_invObliq','var')
            io_utils.WriteSu(filename_invObliq,wavelet_invObliq(:),'dt',dt,'DelayRecordingTime',lagtime*1e3)
        end
    
end