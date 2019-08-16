function clearPathsForSoftRelease(varargin)
% function removes added paths from MATLAB path except for
% $SLIM_APPS/matlab, $SLIM_COMP/matlab, and those that
% match stings specified in the optional cell array of strings
%
% EXAMPLE:
%   clearPathsForSoftRelease
%         clears all except path for clearPathsForSoftRelease
%   clearPathsForSoftRelease({'matlab_etc','matlab_toolboxes'})
%         clears all extras except paths containing 'matlab_etc' and 'matlab_toolboxes'

p = inputParser;
p.addOptional('sfilter',{fileparts(which('clearPathsForSoftRelease'))},@iscell);
p.parse(varargin{:});
sfilter=p.Results.sfilter;

PATH = strsplit(path,':');

    for p=1:length(PATH)
        if ~(find_matlabroot(PATH{p}) ||...
             find_matlabhome(PATH{p}) ||...
             find_SoftRelease(PATH{p}) ||...
             find_other(PATH{p},sfilter)...
            )
            fprintf('removed from path: %s\n',PATH{p})
            rmpath(PATH{p})
        end
    end

end

function test = find_matlabroot(name)
    test = false;
    test = test || ~isempty(strfind(name,matlabroot));
end

function test = find_matlabhome(name)
    test = false;
    filter = fullfile(getenv('HOME'),'matlab');
    test = test || ~isempty(strfind(name,filter));
end

function test = find_SoftRelease(name)
    test = false;
    filter = fullfile(getenv('SLIM_APPS'),'matlab');
    test = test || ~isempty(strfind(name,filter));
    filter = fullfile(getenv('SLIM_COMP'),'matlab');
    test = test || ~isempty(strfind(name,filter));
end

function test = find_other(fname,filter)
    test = false;
    for i=1:length(filter)
        test = test || ~isempty(strfind(fname,filter{i}));
    end
end

