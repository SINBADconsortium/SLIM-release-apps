function ConvertEPS2PDFPNG(dirname)

%% function ConvertEPS2PDFPNG(dirname)
% function to Convert all eps files in the dirname directory to PDF and PNG
% files

%


if nargin < 1
    dirname = './';
end

curdir = pwd;
cd(dirname)

A = dir(['*.eps']);

for i = 1:length(A)
    [path,strfile,ext] = fileparts(A(i).name);
    scommand                = ['unset DYLD_LIBRARY_PATH;convert ' strfile '.eps ' strfile '.png'];
	[status,result]=system(scommand);
    scommand                = ['epstopdf ' strfile '.eps'];
    [status,result]=system(scommand); 
    
end



cd(curdir)