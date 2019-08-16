% This script sets the paths for the data and results. 

curdir = pwd;
basedir = curdir(1:end-48);
datadir = [basedir '/data'];
resultsdir = [basedir '/results/TimeJitteredMarineAcq_OneReceiverGather'];

if ~exist(resultsdir, 'dir')
    mkdir(resultsdir)
end
