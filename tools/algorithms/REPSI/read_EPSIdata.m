function [data dt] = read_EPSIdata(filename, ns, endian)

% Check filetype
if strcmp(filename(end-3:end),'.mat')
    % matlabfile
    filetype = 'matlab';
    disp('Reading from a matlab datafile...')
elseif strcmp(filename(end-3:end),'.bin')
    % pure binary type
    filetype = 'binary';
    if (exist('ns') ~= 1) || isempty(ns)
		error(['Number_of_samples_per_trace must be supplied to read from a binary file!'])
	end
	disp('Reading from a binary datafile...')
elseif strcmp(filename(end-2:end),'.su')
    filetype = 'su';
    disp('Reading from an su datafile...')
else
    error(['Unrecognized filetype for file: ' filename])
end



switch filetype
    
    case 'matlab'
        load(filename,'data');
        checkFor_validDataArray(data)
        
        if size(data, 2) ~= size(data, 3)
            error('n_shot != n_recv in the datafile')
        end
        
        load(filename,'dt');
        
        if ~(exist('dt') == 1)
            dt = [];
        end
        
    case 'binary'
        if (exist('endian') == 1) && ~isempty(endian)
            [data nt ntrace] = io_utils.ReadBinFast(filename, ns, 0, endian);
        else
            [data nt ntrace] = io_utils.ReadBinFast(filename, ns, 0);
        end
        
        try
            checkFor_validDataArray(data)
        catch
            error('found some NaN or Inf values in the data, most likeky cause is wrong endianess')
        end
        
        try
            nrecv = sqrt(ntrace);
            checkFor_validInteger(nrecv)
        catch
            error('square-root of number of traces not an integer, most likely n_shot != n_recv')
        end
        
        nshot = nrecv;
        data = reshape(data, nt, nrecv, nshot);
        dt = [];
        
    case 'su'
        if (exist('endian') == 1) && ~isempty(endian)
            [data nt ntrace dt] = io_utils.ReadSuFast(filename, endian);
        else
            [data nt ntrace dt] = io_utils.ReadSuFast(filename);
        end
        
        try
            checkFor_validDataArray(data)
        catch
            error('found some NaN or Inf values in the data, most likeky cause is wrong endianess')
        end
        
        try
            nrecv = sqrt(ntrace);
            checkFor_validInteger(nrecv)
        catch
            error('square-root of number of traces not an integer, most likely n_shot != n_recv')
        end
        
        nshot = nrecv;
        data = reshape(data, nt, nrecv, nshot);
end