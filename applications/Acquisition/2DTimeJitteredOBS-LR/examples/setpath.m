% This script sets the paths for the data and results. 
% The default is ../../data/ and ../../results.

curdir = pwd;
basedir = curdir(1:end-9);
datadir = [basedir '/data'];
resultsdir = [basedir '/results'];

if ~exist(resultsdir, 'dir')
    mkdir(resultsdir)
end

% add the misc_funcs folder
addpath([basedir '/misc_funcs']);				

