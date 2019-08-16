% This script sets the paths for the data and results. 
% The default is ../../data/ and ../../results.

curdir = pwd;
basedir = curdir(1:end-27);
datadir = [basedir '/data'];
resultsdir = [basedir '/results'];

if ~exist(resultsdir, 'dir')
    mkdir(resultsdir)
end

