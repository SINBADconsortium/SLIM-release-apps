% CLEAN_RSF_HEADER - cleans garbage from RSF header file
%
% CLEAN_RSF_HEADER(header_file_name)
%
function clean_rsf_header(fname)

lines={};
fin=fopen(fname,'r');
while 1
    line=fgetl(fin);
    if ~ischar(line), break, end
    lines{end+1}=line;
end
fclose(fin);

% clear spaces around = sign
lines=strrep(lines,' =','=');
lines=strrep(lines,'= ','=');
% just in case there are two :)
lines=strrep(lines,' =','=');
lines=strrep(lines,'= ','=');

fout=fopen(fname,'w');
for l=1:length(lines)
    fprintf(fout,'%s\n',lines{l});
end
fclose(fout);

end
