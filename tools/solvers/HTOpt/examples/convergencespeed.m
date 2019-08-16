% CONVERGENCESPEED - This script generates an instance of missing tensor interpolation, using seismic data if available, 
% and compares the convergence speed of the available methods in the 
% HT interpolation framework
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%

%% Load seismic data, if available, or generate random data
if exist('slice3d.mat','file') ~= 0    
    load slice3d.mat;
    
    %Low-rank permutation
    fslice = permute(real(fslice),[1 3 2 4]);
    fslice = fslice/norm(vec(fslice));
    dims = size(fslice);
    d = length(dims);
    maxLeafRank = 20;
    maxIntRank = 20;
    dimTree = dimensionTree(dims,maxLeafRank, maxIntRank);
else
    d = 5;  %number of dimensions
    n = 20; %points per dimension
    dims = n*ones(1,d);

    maxLeafRank = 10; 
    maxIntRank = 20; 
    dimTree = dimensionTree(dims,maxLeafRank, maxIntRank);
    x = project(dimTree.randn(),dimTree);
    
    %Normalize data
    [U,B] = dimTree.fromVec(x);
    B{1}{1} = B{1}{1}/norm(vec(B{1}{1}));
    x = dimTree.toVec(U,B);
    
    fslice = dimTree.fullND(x);
end

%% Run experiments
if(exist('saveFigsDir','var')==0)
    error('Need location to save figures');
end
saveFigsDir = [saveFigsDir 'convergencespeed/'];

if(exist(saveFigsDir,'dir')==0)
    mkdir(saveFigsDir);
end
saveDir = saveFigsDir;

opt_methods = {'SD','CG','CG_PRP','GN'};
out_methods = cell(length(opt_methods),1);
for i=1:length(opt_methods)
    rng('default');
    if ~strcmp(opt_methods{i},'GN')
       progTol = 0; % So SD/CG actually runs for sufficiently many iterations 
    else
       progTol = 1e-6; %So GN quits when it starts to stall, near the minimum
    end
    [xsol,rhs,D,out] = HTInterpolation(fslice,dimTree, ...
                                          'subsampling',0.5, ...
                                          'subsample_dims', [1 3], ...
                                          'maxIter',1000, ...
                                          'method',opt_methods{i}, ...
                                          'skipOptimization',false,...
                                          'logFile',[saveDir 'log.txt'],...
                                          'progTol',progTol,...
                                          'tol',1e-6);
    out_methods{i} = out;
end

%% Generate figure

save([saveDir 'convergence.mat'],'out_methods','opt_methods');

fmin = Inf;
maxiter = 0;
for i=1:length(out_methods)
    fmin = min(fmin, min(out_methods{i}.fhist));
    maxiter = max(maxiter, length(out_methods{i}.fhist));
end

saveFig = @(filename) print(gcf, '-depsc',[saveFigsDir filename '.eps']);

figure;     
xdata = 1:maxiter;
ydata = [];
for i=1:length(out_methods)
    ydata = [ydata; [(out_methods{i}.fhist - fmin)' zeros(1,maxiter - length(out_methods{i}.fhist))] ];
end    

colors = {'k','g','b','r'};
handles = [];
for i=1:length(out_methods)
    handles(end+1) = semilogy(xdata,ydata(i,:),colors{i}); hold on; 
end
legend(handles,{'SD','CG', 'CG-PR','GN'});
xlabel('Iteration'); ylabel('Objective value - minimum');
saveFig('convergencespeed');

if exist('metadata','var')
    fid = fopen([saveFigsDir 'log.txt'],'w');
    fprintf(fid,['Script : ' mfilename '\n']);
    fprintf(fid,['Metadata: ' metadata '\n']);
    fclose(fid);
end
