function WriteBin(filename,data)
% Only writes IEEE float32 in default machine endian

segyid = fopen(filename,'w'); 
fwrite(segyid,data,'float32');                 
fclose(segyid);                                  
