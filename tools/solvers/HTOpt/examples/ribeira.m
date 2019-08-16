% RIBEIRA - Script for interpolating the 'Ribeira' dataset. 
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%

rng('default');
load 'data/ribeira.mat';
saveDir = pwd;
if exist(saveDir,'dir')==0
    mkdir(saveDir);
end

saveFigsDir = pwd;
if exist(saveFigsDir,'dir')==0
    mkdir(saveFigsDir);
end

subsampling = [ 0.1, 0.3, 0.5 ];
numTrials = 10;

dims = size(D_ribeira);
tuckRanks1 = [15,15,6];
tuckRanks2 = [65,65,7];

NRMSE = @(X,Xtrue) norm(vec(X) - vec(Xtrue))/((max(vec(Xtrue)) - min(vec(Xtrue))) * sqrt(length(vec(X))));


%%
clear dimTree1 dimTree2;
dimTree1 = dimensionTree(dims,100,268);
dimTree1.rank(3,1,15); dimTree1.rank(3,2,6); dimTree1.rank(2,2,15); dimTree1.rank(2,1,15);

dimTree2 = dimensionTree(dims,100,268);
dimTree2.rank(3,1,65); dimTree2.rank(3,2,7); dimTree2.rank(2,2,65); dimTree2.rank(2,1,65);


%% 
tuck1 = struct;
tuck1.errors = zeros(length(subsampling),numTrials);
tuck1.times = zeros(length(subsampling),numTrials);
tuck2 = struct;
tuck2.errors = zeros(length(subsampling),numTrials);
tuck2.times = zeros(length(subsampling),numTrials);
tuck3 = struct;
tuck3.errors = zeros(length(subsampling),numTrials);
tuck3.times = zeros(length(subsampling),numTrials);
htuck1 = struct;
htuck1.errors = zeros(length(subsampling),numTrials);
htuck1.times = zeros(length(subsampling),numTrials);
htuck2 = struct;
htuck2.errors = zeros(length(subsampling),numTrials);
htuck2.times = zeros(length(subsampling),numTrials);
htuck3 = struct;
htuck3.errors = zeros(length(subsampling),numTrials);
htuck3.times = zeros(length(subsampling),numTrials);

verbosity = 1;

for i=1:length(subsampling)    
    for j=1:numTrials
        trainSet = randperm(prod(dims),round(subsampling(i)*prod(dims)));
        Rtrain = opRestriction(prod(dims),trainSet);
        Rtest = opRestriction(prod(dims),setdiff(1:prod(dims),trainSet));
        
        trainSet = logical(Rtrain' * ones(size(Rtrain,1),1));        
        trainData = D_ribeira;
        trainData(~trainSet) = 0;
        testData = Rtest * vec(D_ribeira);
        
        disp('Tuck 1');
        expTime = tic;
        [U,B] = fitTucker(trainSet, trainData, dims,tuckRanks1,'verbosity',verbosity,'maxIter',200,'progTol',1e-6); 
        expTime = toc(expTime);
        
        X_tuck1 = ttm(B,U,1:3);
        err_t1 = NRMSE(Rtest*vec(ttm(B,U,1:3)),testData);
        tuck1.errors(i,j) = err_t1;
        tuck1.times(i,j) = expTime;
        
        disp('Tuck 2');
        expTime = tic;
        [U,B] = fitTucker(trainSet, trainData, dims,tuckRanks2,'verbosity',verbosity,'maxIter',200,'progTol',1e-6); 
        expTime = toc(expTime);
        
        X_tuck2 = ttm(B,U,1:3);
        err_t2 = NRMSE(Rtest*vec(ttm(B,U,1:3)),testData);
        tuck2.errors(i,j) = err_t2;
        tuck2.times(i,j) = expTime;
        
        disp('Tuck 3');
        numRankSteps = 5;
        expTime = tic;
        [U,B] = fitTuckerRankIncrease(trainSet, trainData, dims, tuckRanks2, numRankSteps, 'verbosity',verbosity,'maxIter',200);
        expTime = toc(expTime);
                
        err_t3 = NRMSE(Rtest*vec(ttm(B,U,1:3)),testData);
        tuck3.errors(i,j) = err_t3;   
        tuck3.times(i,j) = expTime;
        
        X_tuck3 = ttm(B,U,1:3);
        
        disp('Htuck 1');
        expTime = tic;
        x = fitHT(trainSet,trainData,dimTree1,'verbosity',verbosity,'progTol',1e-6,'maxIter',200);
        expTime = toc(expTime);
        
        X_htuck1 = reshape(dimTree1.full(x),dims);
        err_h1 = NRMSE(Rtest*(dimTree1.full(x)),testData);
        htuck1.errors(i,j) = err_h1;
        htuck1.times(i,j) = expTime;
        
        disp('Htuck 2');
        expTime = tic;
        x = fitHT(trainSet,trainData,dimTree2,'verbosity',verbosity,'progTol',1e-6,'maxIter',200);
        expTime = toc(expTime);
        
        X_htuck2 = reshape(dimTree2.full(x),dims);
        err_h2 = NRMSE(Rtest*(dimTree2.full(x)),testData);
        htuck2.errors(i,j) = err_h2;
        htuck2.times(i,j) = expTime;      
        
        disp('Htuck 2');
        expTime = tic;
        x = fitHT(trainSet,trainData,dimTree2,'verbosity',0,'progTol',1e-6,'maxIter',200,'lambda',1e-14);
        expTime = toc(expTime);
        
        X_htuck3 = reshape(dimTree2.full(x),dims);
        err_h3 = NRMSE(Rtest*(dimTree2.full(x)),testData);
        htuck3.errors(i,j) = err_h3;
        htuck3.times(i,j) = expTime;     
        
        if i==1 && j==1
            z = 16;
            slice = @(A) squeeze(A(:,:,z));
            true_slice = slice(D_ribeira);
            imagePlot(true_slice,'cbar',true,'cmap','gray','visible',false);

            saveFigure = @(filename,postfilename) print(gcf,'-depsc',[saveFigsDir filename '-' postfilename '.eps']);
            
            saveFigure('true','');
            
            true_axis = caxis;
            
            interpPlot = @(A) imagePlot(slice(A),'cbar',true,'cmap','gray','visible',false,'coloraxis',true_axis);
            interpPlot(X_tuck2);
            saveFigure('interp-tuck2',['nmsre-' num2str(NRMSE(slice(X_tuck2),slice(D_ribeira)),'%3.3e')]);
            interpPlot(X_tuck3);
            
            saveFigure('interp-tuck3',['nmsre-' num2str(NRMSE(slice(X_tuck3),slice(D_ribeira)),'%3.3e')]);
            
            interpPlot(X_htuck2);
            saveFigure('interp-htuck2',['nmsre-' num2str(NRMSE(slice(X_htuck2),slice(D_ribeira)),'%3.3e')]);                        
            
        end
    end
end

save([saveDir 'ribeira.mat'],'tuck1','tuck2','tuck3','htuck1','htuck2','htuck3');

%%
fid = fopen([saveFigsDir,'ribeira.txt'],'w');
for i=1:length(subsampling)
    fprintf(fid,'------------------------------------\n\n');
    fprintf(fid,['Subsampling ' num2str(subsampling(i)) '\n']);
    fprintf(fid,['Tuck ranks ' num2str(tuckRanks1) '\n']);
    fprintf(fid,['Mean NRMSE ' num2str(mean(tuck1.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD NRMSE ' num2str(std(tuck1.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['Mean Time ' num2str(mean(tuck1.times(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD Time ' num2str(std(tuck1.times(i,:)),'%3.3e') '\n\n']);
    
    fprintf(fid,['Tuck ranks ' num2str(tuckRanks2) '\n']);
    fprintf(fid,['Mean NRMSE ' num2str(mean(tuck2.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD NRMSE ' num2str(std(tuck2.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['Mean Time ' num2str(mean(tuck2.times(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD Time ' num2str(std(tuck2.times(i,:)),'%3.3e') '\n\n']);
    
    fprintf(fid,['Tuck ranks ' num2str(tuckRanks2) ' with continuation\n']);
    fprintf(fid,['Mean NRMSE ' num2str(mean(tuck3.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD NRMSE ' num2str(std(tuck3.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['Mean Time ' num2str(mean(tuck3.times(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD Time ' num2str(std(tuck3.times(i,:)),'%3.3e') '\n\n']);
    
    fprintf(fid,['HTuck ranks ' num2str(tuckRanks1) '\n']);
    fprintf(fid,['Mean NRMSE ' num2str(mean(htuck1.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD NRMSE ' num2str(std(htuck1.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['Mean Time ' num2str(mean(htuck1.times(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD Time ' num2str(std(htuck1.times(i,:)),'%3.3e') '\n\n']);
    
    fprintf(fid,['HTuck ranks ' num2str(tuckRanks2) '\n']);
    fprintf(fid,['Mean NRMSE ' num2str(mean(htuck2.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD NRMSE ' num2str(std(htuck2.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['Mean Time ' num2str(mean(htuck2.times(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD Time ' num2str(std(htuck2.times(i,:)),'%3.3e') '\n\n']);
    
    fprintf(fid,['HTuck ranks ' num2str(tuckRanks2) ' w/ regularization \n']);
    fprintf(fid,['Mean NRMSE ' num2str(mean(htuck3.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD NRMSE ' num2str(std(htuck3.errors(i,:)),'%3.3e') '\n']);
    fprintf(fid,['Mean Time ' num2str(mean(htuck3.times(i,:)),'%3.3e') '\n']);
    fprintf(fid,['STD Time ' num2str(std(htuck3.times(i,:)),'%3.3e') '\n\n']);
end

fclose(fid);

%%
