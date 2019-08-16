function bool = isFile(fileName)
% 
% isFile(fileName)
% 
% Return true if the path points to a file.
% 
import java.io.File;
a=File(fileName);
bool=a.isFile();
end
