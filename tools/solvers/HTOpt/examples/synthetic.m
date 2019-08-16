% SYNTHETICTEST - Generates synthetic HT data + runs HT interpolation for varying
% levels of missing data. 
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%

rng('default');
n = 50; %size of each dimension
d = 4;  %number of dimensions

dims = n * ones(1,d);  %full dimensions
kint = 20; kleaf = 10; %dimension tree ranks

saveFigsDir = pwd;

dimTree = dimensionTree(dims,kleaf,kint); %dimension tree object

x = project(dimTree.randn(),dimTree); %orthogonalized parameters
[U,B] = dimTree.fromVec(x);
B{1}{1} = B{1}{1}/norm(vec(B{1}{1})); 
x = dimTree.toVec(U,B);
X = dimTree.fullND(x); %generate data


% Experiment parameters
numTrials = 2;             
subsamplingRatio = [0.01,0.05,0.1,0.3,0.5];

trainErr = zeros(length(subsamplingRatio),numTrials);
testErr = zeros(length(subsamplingRatio),numTrials);
relErr = zeros(length(subsamplingRatio),numTrials);
solveTime = zeros(length(subsamplingRatio),numTrials);
for i=1:length(subsamplingRatio)
    for j=1:numTrials
        rand_seed = randi(1024);
        [xsol,b,e,out] = HTInterpolation(X,dimTree, ...
                                         'subsampling',subsamplingRatio(i),...
                                         'subsample_dims',1:d,...
                                         'rand_seed',rand_seed, ...
                                         'method','GN',...
                                         'verbosity',1);
        trainErr(i,j) = out.trainErr;
        testErr(i,j) = out.testErr;
        relErr(i,j) = out.relErr;
        solveTime(i,j) = out.solveTime;
    end
end

%%
saveFigure = @(filename) print(gcf,'-depsc',[saveFigsDir filename '.eps']);

mean_errors = mean(relErr,2);
std_errors = std(relErr,1,2)';
figure; 
semilogy(subsamplingRatio, mean_errors','-ok',...
         'LineWidth',3,...
         'MarkerEdgeColor','k',...
         'MarkerFaceColor','g',...
         'MarkerSize',10);

saveFigure('synthetic-meanerr-vs-subsampling');

