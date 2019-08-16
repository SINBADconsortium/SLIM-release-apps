% The example shows on-the-fly extraction of shots/receivers from HT
% parameters. Once you represent your data in compressed HT format either
% in fully sampled scenarios or missing data cases, you can generate the
% shots/receivers from HT parameters without reforming the entire data
% volume, which leads to a substantial reduction in memory costs for
% downstream shot-based algorithms. In this example, we use 3Hz frequency
% slice generated from 3D BG model. 
%
%
% Author: Yiming Zhang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: March, 2018.

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% If you have any questions or issues, pls email
% (yzhang@eoas.ubc.ca)



%% Startup 
close all
clear all
startup 

% Load data
baseDir = pwd;
baseDir = baseDir(1:end-8);
resultsDir = [baseDir 'results/'];
dataDir = [baseDir 'data/'];

rng('default');
if exist(resultsDir,'dir')==0
    mkdir(resultsDir);
end

% If you want to use your own frequency slice, just load it 
% D with size (nrec x X nrec y X nsrc x X nsrc y)
load([dataDir '/BG_3Hz.mat']);

%% In the fully sampled case, you can truncate your data D in compressed HT parameters

% Permute your data into noncanical organization (receiver x, source x, receiver y, source y)
D1               = permute(D,[1 3 2 4]);

% Truncates to HT format, given the relative error bound 1e-3, and you can change the error
% tolerance as you want 
[dimTree1,x1]     = truncate_ht(D1,1e-3,1,0);

% A shot generation from vecterized HT parameters
geometry.so = [0, 0];
geometry.ro = [0, 0];
geometry.sd = [200, 200];
geometry.rd = [50, 50];
geometry.sn = [49, 49];
geometry.rn = [196, 196];
index1 = [55];
mode1  = '1';

shot1 = slicing(dimTree1, x1, index1, geometry, mode1);

% plot 
figure;imagesc(reshape(real(shot1),196,196));colormap seiscol
figure;imagesc(real(squeeze(D(:,:,6,2)))); colormap seiscol


%% In the missing entries case, you can interpolate HT format 

% src rank
ksrc = 49;

% (recx srcx) rank
kint = 470;

% rec rank
krec = 55;

% Number of sources and receivers along one direction  
nrecs = 196; nsrcs = 49;

% Percentage sources to throw away
subPercent = 0.75;

% Maximum iterations
maxIter = 30;

% Available methods : 'SD','CG','CG_PR','GN'
method = 'GN';

% 0 - silent, 1 - iteration, 2 - verbose
verbosity = 1;

% Remove [recx recy] pairs
subDims =[1 2];

% Regularization parameter - mostly useful when subPercent >= 0.9
lambda = 0;

% Dimensions of D
dims = [nrecs,nrecs,nsrcs,nsrcs];

% Permuted dimensions
% P : (recx, srcx, recy,srcy) -> (recx, recy,srcx, srcy)
% P': (recx, recy,srcx, srcy)-> (recx, srcx, recy,srcy)
pdims = [nrecs,nsrcs, nrecs,nsrcs];
P = opPermute(pdims, [1 3 2 4]);

parallelMode= false;
% Generate subsampling operators
[Rtrain,Rtest,Jtrain,Jtest] = ndimSubsampling(dims,subDims,subPercent,0,parallelMode);

% Subsampling operator
A = Rtrain * P;

% For point-removal sampling, it's faster to store 1s where there
% are data points (vector e)

e = A' * ones(size(A,1),1);
e = logical(e);

% Subsample data points, in the [rec x, src x, rec y, src y] domain
bsub = P' * vec(D);
bsub(~e) = 0;


%Noncanonical dimension tree
dimTree2 = dimensionTree(pdims,max(ksrc,krec),kint);
dimTree2.rank(3,1,krec); dimTree2.rank(3,2,ksrc);
dimTree2.rank(3,3,krec); dimTree2.rank(3,4,ksrc);

% Initial point
x0 = project(dimTree2.randn(),dimTree2);
[U,B] = dimTree2.fromVec(x0);
B{1}{1} = B{1}{1}/norm(vec(B{1}{1}));
x0 = dimTree2.toVec(U,B);


% Solve interpolation problem
funObj = @(x) LSLinearMisfit(x, bsub, ~e);

x2 = minFunc_hTuck(funObj, dimTree2, ...
                  'maxIter',maxIter,...
                  'method',method,...
                  'x0',x0,...
                  'lambda',lambda,...
                  'verbosity',verbosity,...
                  'distributed',parallelMode);
              
% A shot generation from vecterized HT parameters     
index2 = [55];
mode2  = '1';

shot2 = slicing(dimTree2, x2, index2, geometry, mode2);

% plot 
figure;imagesc(reshape(real(shot2),196,196));colormap seiscol
figure;imagesc(real(squeeze(D(:,:,6,2)))); colormap seiscol
              
 %% Follow the idea of extraction of shots/receivers from compressed HT parameter, you can compute D*w or D^H*w if you provide the probing vector

 % Probing vectors
 v1 = randn(nsrcs*nsrcs,1);
 v2 = randn(nrecs*nrecs,1);
 
 % Form the conventional data matrix 
 D2 = reshape(D, [nrecs*nrecs, nsrcs*nsrcs]);
 
 % Compute D*v1 and D^H*v2 using D
 d1 = D2*v1;
 d2 = D2'*v2;
 
 % Using HT parameters
 mode = ['1','2'];
 d3   = Dv_HT(x1,dimTree1,v1,mode(2));
 d4   = Dv_HT(x1,dimTree1,v2,mode(1));

 % plot
 figure;imagesc(reshape(real(d1), nrecs, nrecs)); colormap seiscol
 figure;imagesc(reshape(real(d3), nrecs, nrecs)); colormap seiscol
 figure;imagesc(reshape(real(d2), nsrcs, nsrcs)); colormap seiscol
 figure;imagesc(reshape(real(d4), nsrcs, nsrcs)); colormap seiscol
 
 %% Save results
 save([resultsDir 'results.mat'],'dimTree1','x1','dimTree2','x2','shot1','shot2','d1','d2','d3','d4');
 
