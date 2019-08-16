% INTERP4D - Interpolates 4D seismic data with missing receivers
% using the Hierarchical Tucker format. Uses dense linear algebra
% routines that form the entire intermediate tensor at each
% iteration. Serial and parallel routines provided.
%
% Author: Curt Da Silva (curtd@math.ubc.ca)
% Date: March 2014

%% Load data
baseDir = pwd;
baseDir = baseDir(1:end-8);
resultsDir = [baseDir 'results/'];
dataDir = [baseDir 'data/'];

rng('default');
if exist(resultsDir,'dir')==0
    mkdir(resultsDir);
end

% If you want to use your own frequency slice, just load it in matlab and call
% it D (with size nsrcs x nsrcs x nrecs x nrecs)
load([dataDir '/slice3d.mat']);

% We organize the data as nrecs x nrecs x nsrcs x nsrcs, because
% matlab parallelizes its data over the outer dimensions and we can
% parallelize over shots easier this way.
D = D/norm(vec(D)); D = permute(D,[3 4 1 2]);

nsrcs = size(D,1); nrecs = size(D,3);

% Toggle parallel distribution or not
% Will only really be useful when D is too big to fit in memory of
% a single machine, since dense linear algebra routines are multithreaded already
parallelMode = false;

if parallelMode
    parpool_close()
    parpool_open(4);
    D = matricize(D,[1 2]);
    spmd,
        D = codistributed(D);
    end
end

%% Experiment parameters
%HT parameters
    
% src rank
ksrc = 20; 
% (srcx recx) rank
kint = 20;    
% rec rank
krec = 20;    

% Percentage receivers to throw away
subPercent = 0.75;                

% Maximum iterations
maxIter = 20;       

% Available methods : 'SD','CG','CG_PR','GN', see minFunc_hTuck.m for more options
method = 'GN';

% 0 - silent, 1 - iteration, 2 - verbose
verbosity = 1;

% Remove [recx recy] pairs
subDims = [1 2];

% Regularization parameter - mostly useful when subPercent >= 0.9
lambda = 0;

%% Subsampling
% Dimensions of D
dims = [nrecs,nrecs,nsrcs,nsrcs];

% Permuted dimensions
% P : (recx, srcx, recy,srcy) -> (recx, recy,srcx, srcy)
% P': (recx, recy,srcx, srcy)-> (recx, srcx, recy,srcy)
pdims = [nrecs,nsrcs, nrecs,nsrcs];    
P = opPermute(pdims, [1 3 2 4]);

% Generate subsampling operators
[Rtrain,Rtest,Jtrain,Jtest] = ndimSubsampling(dims,subDims,subPercent,0,parallelMode);

% Subsampling operator
A = Rtrain * P;

% For point-removal sampling, it's faster to store 1s where there
% are data points (vector e)
if parallelMode   
    
    btmp = Rtrain * pSPOT.utils.distVectorize(D);    
    spmd
        f = codistributed.ones(length(btmp),1,getCodistributor(btmp));
    end
    clear btmp;
    e = A' * f;
    clear f;
    e = logical(e);
    bsub = P' * pSPOT.utils.distVectorize(D);
    bsub(~e) = 0;
else   
    e = A' * ones(size(A,1),1);
    e = logical(e);
    
    % Subsample data points, in the [src x, rec x, src y, rec y] domain
    bsub = P' * vec(D);
    bsub(~e) = 0;
end

%% Create dimension tree
%Noncanonical dimension tree
dimTree = dimensionTree(pdims,max(ksrc,krec),kint);
dimTree.rank(3,1,krec); dimTree.rank(3,2,ksrc);
dimTree.rank(3,3,krec); dimTree.rank(3,4,ksrc);

% Initial point
x0 = project(dimTree.randn(),dimTree);
[U,B] = dimTree.fromVec(x0);
B{1}{1} = B{1}{1}/norm(vec(B{1}{1}));
x0 = dimTree.toVec(U,B);

if parallelMode
    %Need to make sure that bsub, e are distributed identically to
    %the distributed HT tensor. 
    Xfull = dimTree.fullDist(x0);
    spmd
        e = redistribute(e,getCodistributor(Xfull));
        bsub = redistribute(bsub,getCodistributor(Xfull));
    end
    clear Xfull; 
end

% Solve interpolation problem
funObj = @(x) LSLinearMisfit(x, bsub, ~e);    

x = minFunc_hTuck(funObj, dimTree, ...
                  'maxIter',maxIter,...
                  'method',method,...
                  'x0',x0,...
                  'lambda',lambda,...
                  'verbosity',verbosity,...
                  'distributed',parallelMode);

% Arrange everything in to the recx x recy x srcx x srcy domain,
% for plotting purposes.
if parallelMode
    X = dimTree.fullDist(x);
    D = reshape(gather(D),dims);
    X = P' * gather(X);
    bsub = reshape(P*gather(bsub),dims);
    e_train = P * gather(e);
else    
    X = reshape(P * dimTree.full(x),dims);
    e_train = logical(Rtrain' * ones(size(Rtrain,1),1));
    bsub = reshape(P*bsub,dims);
end

% Compute recovery 
test_D = vec(D); test_D(e_train) = 0;
test_X = vec(X); test_X(e_train) = 0;
test_snr = SNR(test_D,test_X);
overall_snr = SNR(vec(D),vec(X));


%save([resultsDir 'results.mat'],'D','X','bsub','test_snr','overall_snr');

