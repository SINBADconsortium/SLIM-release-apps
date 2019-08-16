% INTERP4D - Interpolates 4D seismic data with missing receivers
% using the Hierarchical Tucker format. Uses sparse linear algebra to implicitly form the full tensor at each iteration. More efficient than interp4D.m only when the tensor size is much larger than the amount of data available.
% 
% - Serial and parallel routines provided, but parallel is recommended
% - Can only handle real data for the time being 
%
% Author: Curt Da Silva (curtd@math.ubc.ca)
% Date: March 2014

%% Load data
baseDir = pwd;
baseDir = baseDir(1:end-8);
resultsDir = [baseDir 'results/'];
dataDir = [baseDir 'data/'];

if exist(resultsDir,'dir')==0
    mkdir(resultsDir);
end

rng('default');

% If you want to use your own frequency slice, just load it in matlab and call
% it D (with size nsrcs x nsrcs x nrecs x nrecs)
load([dataDir 'slice3d.mat']);

% We organize the data as nrecs x nrecs x nsrcs x nsrcs, because
% matlab parallelizes its data over the outer dimensions and we can
% parallelize over shots easier this way.
D = real(D); D = D/norm(vec(D)); D = permute(D,[3 4 1 2]);

nsrcs = size(D,1); nrecs = size(D,3);

% Toggle parallel distribution or not
% Will only really be useful when D is too big to fit in memory,
% since dense linear algebra routines are multithreaded already
parallelMode = true;

if parallelMode
    assert(parpool_size==8,'size of parallel pool should be 8');
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
maxIter = 50;       

% Available methods : 'SD','CG','CG_PR','GN', see minFunc_hTuck.m for more options
method = 'GN';

% 0 - silent, 1 - iteration, 2 - verbose
verbosity = 1;

% Remove [recx recy] pairs
subDims = [1 2];

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

% I = d x numdatapts matrix of d-dimensional indices where there is data
if parallelMode   
    bsub = P' * pSPOT.utils.distVectorize(D);
    btmp = Rtrain * pSPOT.utils.distVectorize(D);    
    spmd
        f = codistributed.ones(length(btmp),1,getCodistributor(btmp));
    end
    clear btmp;
    e = A' * f;
    clear f;
    e = logical(e);
    I = find(e);
    bsub_sparse = bsub(I);  
    
    % Construct I in a distributed fashion
    spmd
        Iloc = getLocalPart(I);        
        Ifullloc = uint32(idx1d2nd(pdims,Iloc));        
        codist = getCodistributor(I);
        partition = codist.Partition;
        codist = codistributor1d(2,partition,[length(pdims),sum(partition)]);
        I = codistributed.build(Ifullloc,codist); 
    end
    
    % I, bsub should have the same distribution
    spmd
        part = getCodistributor(I);
        part = part.Partition;
        codist = codistributor1d(1,part,[sum(part),1]);
        bsub_sparse = redistribute(bsub_sparse,codist);
    end   
    
else   
    e = A' * ones(size(A,1),1);
    e = logical(e);
    I = uint32(idx1d2nd(pdims,find(e)));
    
    % Subsample data points, in the [src x, rec x, src y, rec y] domain
    bsub = P' * vec(D);
    bsub(~e) = 0;    
    bsub_sparse = bsub(e); 
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

% Solve interpolation problem
funObj = @(x) LSMisfitHTmex(x, I,bsub_sparse,dimTree);   
x = minFunc_hTuck(funObj, dimTree, ...
                  'maxIter',maxIter,...
                  'method',method,...
                  'x0',x0,...
                  'verbosity',verbosity,...
                  'dimtreeFunc',true);


% Arrange everything in to the recx x recy x srcx x srcy domain,
% for plotting purposes.
if parallelMode
    X = P * dimTree.full(x);
    D = reshape(gather(D),dims);    
    e_train = P * gather(e);
    bsub = reshape(P * gather(bsub),dims);
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


save([resultsDir 'results.mat'],'X','bsub','test_snr','overall_snr');



