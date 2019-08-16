function [a,o,d,n] = odnread(fname)
% read regularly gridded data from .rsf-like file. see also odnwrite.
%
% use:
%   [a,o,d,n] = odnread(fname)
%
% input:
%   fname - name of file
%
% output:
%   a       - data
%   {o,d,n} - grid parameters
%
% odn can be converted to grid with odn2grid utility
%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

try 
    fidh = fopen(fname,'r');
    fidb = fopen([fname '@'],'r');
catch
    fprintf(2,'File not found\n');
end

% read header information
k = 1;
line  = fgetl(fidh); 
tmp = sscanf(line,['n' num2str(k) '=%d d' num2str(k) '=%f o' num2str(k) '=%f']);
while ~isempty(tmp)
    if ~isempty(tmp)
        n(k) = tmp(1); d(k) = tmp(2); o(k) = tmp(3);
        k = k + 1;
    end
    line  = fgetl(fidh); 
    tmp = sscanf(line,['n' num2str(k) '=%d d' num2str(k) '=%f o' num2str(k) '=%f']); 
end
line = fgetl(fidh);
type = sscanf(line,'data_type=%s');

% read binary part
a    = fread(fidb,'float');

% check for imaginary part 
if strcmp(type,'complex')
    a = a(1:length(a)/2) + 1i*a(length(a)/2+1:end);
end

% close files
fclose(fidh);fclose(fidb);