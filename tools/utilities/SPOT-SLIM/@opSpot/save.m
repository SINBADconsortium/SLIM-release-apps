function save(obj,fileName,overwrite)
%SAVE Summary of this function goes here
%   Detailed explanation goes here

if(nargin == 2)
    overwrite = 0;
end

if (overwrite == 1)
    if (exist(fileName))
        delete(fileName);
    end
else
    assert(~exist(fileName),'file already exists')
end

save(fileName,'obj')
end

