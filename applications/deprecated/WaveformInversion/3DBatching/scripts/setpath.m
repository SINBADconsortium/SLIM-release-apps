% This script sets the paths for the results and data. 
% The default is ../results and ../data/.

curdir     = pwd;
basedir    = curdir(1:end-8);

datadir    = [basedir '/data'];
resultsdir = [basedir '/results'];

if ~exist(resultsdir,'dir')
    mkdir(resultsdir)
end

setenv('OMP_NUM_THREADS', '10');
