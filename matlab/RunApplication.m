function RunApplication(AppDir,AppStart,ExpDir,ExpStart,Applications)
% SINBAD application helper
%
% RunApplication runs the application in the SINBAD software release
% ensuring that the proper MATLAB path environment is loaded in both
% interactive and batch job (see batch help for running functions
% with arguments).
%
% RunApplication(AppDir,AppStart,ExpDir,ExpStart,Applications)
% where:
%   AppDir is the string with path to main application directory
%       relative to $SLIM_APPS_RUNS
%   AppStart is a logical flag if to run startup.m from AppDir
%   ExpDir is the string with path to experiment subdirectory
%       relative to AppDir
%   ExpStart is a logical flag if to run startup.m from ExpDir
%   Applications is a cell array of strings with MATLAB script names to run in order
%
% Notes:
%   - this function will remove obsolete paths from MATLAB path
%   - parallel pool is not created by RunApplication and must be
%     open prior to using RunApplications
%
% Examples:
%   running two scripts where the later depends on the results from former
%     RunApplication('applications/Processing/SparsityPromotingDenoising',...
%       true,'examples',false,{'getCurveletCoeff','denoiseSeismicData'})
%   running the script if experiment sub-directory is the same as main directory
%     RunApplication('applications/Imaging/WRimaging',...
%       true,'./',false,{'WRImaging_example'})
%   running aplication if startup.m is in experiment sub-directory
%     RunApplication('applications/Modeling/TTITimeModeling',...
%       false,'scripts',true,{'Modelling_TTI_2D'})
%

    clearPathsForSoftRelease
    olddir=pwd;
    home_dir=getenv('HOME');
    slim_apps_runs=getenv('SLIM_APPS_RUNS');
    assert(~isempty(slim_apps_runs),'SLIM_APPS_RUNS environment is not set. Check environment.sh source')
    app_dir=fullfile(slim_apps_runs,AppDir);
    exp_dir=fullfile(slim_apps_runs,AppDir,ExpDir);

    changedir(app_dir,AppStart)
    changedir(exp_dir,ExpStart)

    for a=1:length(Applications)
        app = Applications{a};
        fprintf('----- starting %s\n',app);
        disp(app);
        eval(app)
        fprintf('----- done %s\n',app);
    end

    cd(olddir);
    fprintf('----- all done\n');
end

function changedir(MyDir,MyStart)
    assert(isdir(MyDir),'Directory %s does not exist',MyDir);
    fprintf('----- Changing directory to %s\n',MyDir);
    cd(MyDir)
    load_startup(MyStart)
end

function load_startup(flag)
    if flag
        assert(isFile('startup.m'),'startup.m not found in current directory');
        if exist('startup.m','file')
            fprintf('----- Loading startup in %s\n',pwd);
            startup
        end
    end
end
