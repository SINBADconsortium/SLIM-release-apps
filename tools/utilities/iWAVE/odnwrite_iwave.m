function [] = odnwrite_iwave(fname,a,o,d,n,type)
% write regularly gridded data to .rsf-like file. see also odnread.
%
% use:
%   [] = odnwrite(fname,a,o,d,n)
%
% input:
%   fname   - filename
%   a       - data
%   {o,d,n} - regular grid parameters: z = o(1) + [0:n(1)-1]*d(1), etc.
%
% [o,d,n] can be produced by grid2odn(z,x,y) utility
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

if nargin<6
    if isreal(a)
        type = 'real';
    else
        type = 'complex';
    end
end

try
    fidh = fopen(fname,'w');
    fidb = fopen([fname '@'],'w');
catch
    fprintf(2,'Could not open file\n');
end

% vectorize data
a = a(:);

% check for imaginary part, write binary part
if isreal(a)
    fwrite(fidb,a,'float');
else
    fwrite(fidb,[real(a) imag(a)],'float');
end

% write header
for k = 1:length(n)
    fprintf(fidh,'n%d=%d d%d=%f o%d=%f\n',k,n(k),k,d(k),k,o(k));
end
fprintf(fidh,'data_format=native_float\n');
fprintf(fidh,'data_type=%s\n',type);
fprintf(fidh,'in=%s\n',[fname '@']);

% close files
fclose(fidh);
fclose(fidb);
