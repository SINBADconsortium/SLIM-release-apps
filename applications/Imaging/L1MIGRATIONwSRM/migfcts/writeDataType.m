function writeDataType(filename,data,variable_name)
% Syntax:
% writeDataType(filename,data,variable_name)
%
% Description:
% write 'data' to file 'filename'. If you are saving a matlab file, the variable
% name can also be defined with the optional input "variable_name". It currently
% supports SU/RSF/binary/MAT files given ".su/.rsf/.bin/.mat" suffix in file
% name.
%
% Input list:
% filename: name of file you want to save your results to
% data: the result you want to write
% variable_name: for matlab files only. Designate variable name for matlab
% 	variables. Default is "data".
%
% Output list: None. Save data to files.
%
% Author: Ning Tu (Major contribution from Tim Lin)
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if strcmp(filename(end-3:end),'.mat')
    % matlabfile
    filetype = 'matlab';
    disp('Writing to a matlab datafile...')
elseif strcmp(filename(end-3:end),'.bin')
    % pure binary type
    filetype = 'binary';
	disp('Writing to a binary datafile...')
elseif strcmp(filename(end-3:end),'.rsf')
    filetype = 'rsf';
    disp('Writing to a rsf datafile...')
elseif strcmp(filename(end-2:end),'.su')
    filetype = 'su';
    disp('Writing to a su datafile...')
else
    error(['Unrecognized filetype for output file: ' filename])
end

switch filetype
    
    case 'matlab'
        if nargin < 3
            save(filename,'data');
        else
            eval([variable_name ' = data;'])
            eval(['save(''' filename ''',''' variable_name ''')'])
        end
        
    case 'binary'
        io_utils.WriteBin(filename,data)
        
    case 'rsf'
        rsf_write_all(filename,{},data)
        
    case 'su'
        dim1 = size(data,1);
        dim2 = size(data,2);
        dim3 = size(data,3);
        data = reshape(data, dim1, dim2*dim3);
        
        % Make an array for the correct FieldRecord value for the shot records
        groupId = zeros(dim2*dim3,1);
        for k = 1:dim3
            groupId((k-1)*dim2+1:k*dim2) = k;
        end
        
        traceId = zeros(dim2*dim3,1);
        for k = 1:dim3
            traceId((k-1)*dim2+1:k*dim2) = [1:dim2];
        end
        
        io_utils.WriteSu(filename,data,'FieldRecord',groupId.','TraceNumber',traceId.')
end