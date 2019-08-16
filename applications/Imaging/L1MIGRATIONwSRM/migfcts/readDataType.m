function data = readDataType(filename, variablename, endian, ns)
% Syntax:
% data = readDataType(filename, variablename, endian, ns)
% 
% Description:
% load data from file. It currently supports SU/RSF/binary/MAT files given 
% ".su/.rsf/.bin/.mat" suffix in file name.
% 
% Input list:
% filename: name of file you read data from
% variablename: necessary if file is matlab file, otherwise set empty
% endian: optional, can help read binary and SU files
% ns: number of samples: necessary if reading binary files
%
% Output list:
% data: the data you read
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

% Check filetype
if strcmp(filename(end-3:end),'.mat')
    % matlabfile
    filetype = 'matlab';
    disp('Reading from a matlab datafile...')
elseif strcmp(filename(end-3:end),'.bin')
    % pure binary type
    filetype = 'binary';
    if isempty(ns)
        error(['Number_of_samples_per_trace must be supplied to read from a binary file!'])
    end
    disp('Reading from a binary datafile...')
elseif strcmp(filename(end-2:end),'.su')
    filetype = 'su';
    disp('Reading from an su datafile...')
elseif strcmp(filename(end-3:end),'.rsf')
    filetype = 'rsf';
    disp('Reading from RSF datafile...')
    if ~exist('rsf_dim')
        error('RSF matlab API has not been installed correctly')
    end
else
    error(['Unrecognized filetype for file:' filename])
end

% check input args
if not(exist('endian','var'))
    endian = [];
end
if not(exist('ns','var'))
    ns = [];
end

% read files
switch filetype
    
  case 'matlab'
    load(filename,variablename);
    if not(isequal(variablename,'data'))
        eval(['data = ' variablename ';']);
        eval(['clear ' variablename]);
    end
    checkFor_validDataArray(data)
  case 'binary'
    if isempty(ns)
        error(['Number of samples in the 1st dimension must be designated for binary files']);
    end
    if ~isempty(endian)
        [data nt ntrace] = io_utils.ReadBinFast(filename, ns, 0, endian);
    else
        [data nt ntrace] = io_utils.ReadBinFast(filename, ns, 0);
    end
    
    try
        checkFor_validDataArray(data)
    catch
        error('found some NaN or Inf values in the data, most likely cause is wrong endianess')
    end
    
    data = reshape(data,nt,ntrace);
  case 'su'
    if ~isempty(endian)
        [data nt ntrace] = io_utils.ReadSuFast(filename, endian);
    else
        [data nt ntrace] = io_utils.ReadSuFast(filename);
    end
    
    try
        checkFor_validDataArray(data)
    catch
        error('found some NaN or Inf values in the data, most likely cause is wrong endianess')
    end
    
    data = reshape(data,nt,ntrace);
  case 'rsf'
    data = rsf_read_all(filename);
    
    checkFor_validDataArray(data)
end