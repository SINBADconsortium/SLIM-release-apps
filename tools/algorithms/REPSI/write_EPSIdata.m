function write_EPSIdata(filename,primary,dt)
% Writes datacubes

% Check filetype
if strcmp(filename(end-3:end),'.mat')
    % matlabfile
    filetype = 'matlab';
    disp('Writing primaries to a matlab datafile...')
elseif strcmp(filename(end-3:end),'.bin')
    % pure binary type
    filetype = 'binary';
	disp('Writing primaries to a binary datafile...')
elseif strcmp(filename(end-2:end),'.su')
    filetype = 'su';
    disp('Writing primaries to an su datafile... (WARNING: could take several minutes)')
else
    error(['Unrecognized filetype for output file: ' filename])
end



switch filetype
    
    case 'matlab'
        save(filename,'primary','dt')
        
    case 'binary'
        nt = size(primary,1);
        nr = size(primary,2);
        ns = size(primary,3);
        primary = reshape(primary, nt, nr*ns);
        
        io_utils.WriteBin(filename,primary)
        
    case 'su'
        nt = size(primary,1);
        nr = size(primary,2);
        ns = size(primary,3);
        primary = reshape(primary, nt, nr*ns);
        
        % Make an array for the correct FieldRecord value for the shot records
        groupId = zeros(nr*ns,1);
        for k = 1:ns
            groupId((k-1)*nr+1:k*nr) = k;
        end
        
        traceId = zeros(nr*ns,1);
        for k = 1:ns
            traceId((k-1)*nr+1:k*nr) = [1:nr];
        end
        
        io_utils.WriteSu(filename,primary,'dt',dt,'FieldRecord',groupId.','TraceNumber',traceId.')
    
end