% This script sets the paths for the data and results. 
% The default is ../data/ and ../results.

curdir     = pwd;
basedir    = curdir(1:end-20);
datainterp = [basedir '/data/Interpolation/'];
datareg    = [basedir '/data/Regularization/'];
resultsdir = [basedir '/results/'];

if ~exist(resultsdir, 'dir')
    mkdir(resultsdir)
end