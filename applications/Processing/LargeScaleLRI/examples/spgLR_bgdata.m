%% Print figure parameters

srcx = 30; srcy = 30;

% Verbosity level
verbosity = 1;

% How often to display output
outputFreq = 50;

% Which experiment to run, see spglr_experiments/ directory for details
experiment = 1;

% Random seed
randseed = 1752834;

%% Load data
baseDir = pwd;
baseDir = baseDir(1:end-8);
resultsDir = [baseDir 'results/'];
dataDir = [baseDir 'data/'];
rng('default');
if exist(resultsDir,'dir')==0
    mkdir(resultsDir);
end

assert(exist('experiment','var')==1,'Need experiment variable');
assert(parpool_size()>0,'Need open parallel pool');

numWorkers = parpool_size();

% Get experiment parameters
eval(['spgLR_experiment' num2str(experiment)]);    
v2struct(params);

rng(rand_seed);

% Load data
switch freqIndex
  case 75
    dload([dataDir 'bgdata_freq75.mat']);
  case 125
    dload([dataDir 'bgdata_freq125.mat']);
end

% the original data b has dimensions [nrecs,nrecs,nsrcs,nsrcs]
%
% the provided b has been matricized, so b = matricize(b_original,[1 3],[2 4]);
%
nrecs = dims(1); nsrcs = dims(3);

m = nsrcs*nrecs; n = nsrcs*nrecs;

% Set up results directories
saveDir = resultsDir;

%% Set up permutation, subsampling operators
pdims = [nrecs,nsrcs,nrecs,nsrcs];
P = oppPermute(pdims,[1 3 2 4]);
J = sort(randperm(nrecs^2,round((1-subPercent)*nrecs^2)),'ascend');

spmd,    
    codist = codistributor1d(2,codistributor1d.unsetPartition,[nrecs^2,nsrcs^2]);
    e = false([nrecs^2,nsrcs^2],codist);
    eloc = getLocalPart(e);
    eloc(J,:) = true;  
    e = codistributed.build(eloc,getCodistributor(e));
end
e = pSPOT.utils.distVectorize(e);
e = pSPOT.utils.distVec2distArray( P' * e, [m,n]);
e = logical(e); 

%% b is already in permuted form, so just subsample it
spmd,
    e = redistribute(e,getCodistributor(b));
    bloc = getLocalPart(b); eloc = getLocalPart(e);
    bloc(~eloc) = 0;
    b = codistributed.build(bloc,getCodistributor(b));
    meanloc = sum(abs(full(vec(bloc)))); meanMag = pSPOT.utils.global_sum(meanloc); 
    bloc=[]; eloc = [];
end

% Mean magnitude of each entry
meanMag = meanMag{1}/((1-subPercent)*numel(b));

%% Actual optimization code
bNorm = norm(b,'fro');
opts = struct; 
opts.maxIter = maxIter; 
opts.scale = meanMag;
opts.verbosity = verbosity; 
opts.outputFreq = outputFreq; 
if exist('Lest','var')
    opts.Linit = Lest; opts.Rinit = Rest;
end

tic;
[Lest,Rest,res,spgLRTime] = spgLR(b,e,rank,sigma*bNorm,opts);
spgLRtime = toc;

%% Save results
saveFile = [saveDir 'spgLR_experiment' num2str(experiment) '.mat'];
dsave(saveFile,'Lest','Rest','spgLRTime','e');

