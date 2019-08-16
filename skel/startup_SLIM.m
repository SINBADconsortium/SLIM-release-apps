% Software Release addons
if ~isempty(getenv('SLIM_COMP'))
    slim_comp=fullfile(getenv('SLIM_COMP'),'matlab');
    if isdir(slim_comp);
        addpath(slim_comp)
        fprintf('Added SLIM_COMP/matlab to path\n')
    end
else
    fprintf('SLIM_COMP undefined. Skipped adding its matlab path.\n')
end
if ~isempty(getenv('SLIM_APPS'))
    slim_apps=fullfile(getenv('SLIM_APPS'),'matlab');
    if isdir(slim_apps)
        addpath(slim_apps);
        fprintf('Added SLIM_APPS/matlab to path\n')
    end
else
    error('Software:Release:Error','FATAL ERROR: missing environment SLIM_APPS');
end
