function bool = isDirectory(dirName)
% 
% isDirectory(dirName)
% 
% Return true if the path points to a directory.
% 
import java.io.File;
a=File(dirName);
bool=a.isDirectory();
end
