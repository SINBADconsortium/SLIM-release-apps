function flag = addSRpath(srhome,spath,rcrsv)
% add path to MATLAB path, or move it to the front if it was added before
% addSRpath(srhome,spath,rcrsv)
% where:
%   srhome is string 'SLIM_APPS' or 'SLIM_COMP'
%   spath is path to add relative to above
%   rcrcs if to use recursive genpath
%
    flag = true;
    basepath = getenv(srhome);
    if isempty(basepath) or ~isdir(basepath)
	    fprintf('FATAL ERROR: %s environment or its path does not exist\n',srhome);
	    fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
        flag=false;
	    return;
    end
    mypath = fullfile(basepath,spath);
    assert(isdir(mypath'),'FATAL ERROR: %s directory does not exist in %s',spath,basepath);

    myrmpath(mypath, rcrsv)
    myaddpath(mypath, rcrsv)

end

function myrmpath(mydir, gflag)
    PATH = strsplit(path,':');
    if gflag
        mydirs = strsplit(genpath(mydir),':');
        for p=1:length(PATH)
            for m=1:length(mydirs)
                if strcmp(PATH{p},mydirs{m})
                    rmpath(mydirs{m});
                end
            end
        end
    else
        for p=length(PATH):-1:1
            if strcmp(PATH{p},mydir)
                rmpath(PATH{p});
            end
        end
    end
end

function myaddpath(mydir, gflag)
    if gflag
        mydirs = strsplit(genpath(mydir),':');
        for m=length(mydirs):-1:1
            addpath(mydirs{m});
        end
    else
        addpath(mydir);
    end
end
