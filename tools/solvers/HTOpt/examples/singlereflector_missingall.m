% SINGLEREFLECTOR_MISSINGRECS - Interpolates seismic data with
% missing receivers. Compares the Tucker approach vs the Hierarchical
% Tucker approach.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca

%% Load data
load '~/code/HTpaper/slice3d.mat';
rng('default');

% Location to save figures/tables
if exist('saveFigsDir','var')==0
    error('Need location to save figures');
else
    saveFigsDir = [saveFigsDir 'singlereflector-missingall/'];
end

saveDir = saveFigsDir;

if exist(saveFigsDir,'dir')==0
    mkdir(saveFigsDir);
end

%% Experimental parameters
subsampling = [ 0.1, 0.3, 0.5 ];
numTrials = 5;

% (Src x, Rec x) matricization/dimension tree
fslice = real(fslice); fslice = fslice/norm(vec(fslice));
fslice = permute(fslice,[1 3 2 4]);

%Tucker parameters
dims = size(fslice);
tuckRanks1 = 20*ones(1,4);
tuckRanks2 = 30*ones(1,4);

range = @(y) max(vec(y)) - min(vec(y));

NRMSE = @(X,Xtrue) norm(vec(X) - vec(Xtrue))/( ( range(real(Xtrue)) + range(imag(Xtrue)) )  * sqrt(numel(X)) );
% HTucker parameters
clear dimTree1 dimTree2 dimTree3;
kint1 = 60; kleaf1 = 20;
dimTree1 = dimensionTree(dims,kleaf1,kint1);
kint2 = 80; kleaf2 = 20;
dimTree2 = dimensionTree(dims,kleaf2,kint2);
kint3 = 80; kleaf3 = 30;
dimTree3 = dimensionTree(dims,kleaf3,kint3);

%% Experiment results
% Tuck1 - tucker format with ranks tuckRanks1
% Tuck2 - tucker format with ranks tuckRanks2
% Tuck3 - tucker format with ranks tuckRanks2 + rank continuation
% Htuck1 - HT format with ranks specified by dimTree1
% Htuck2 - HT format with ranks specified by dimTree2
% Htuck3 - HT format with ranks specified by dimTree3
tuck1 = struct;
tuck1.snr = zeros(length(subsampling),numTrials);
tuck1.times = zeros(length(subsampling),numTrials);
tuck2 = struct;
tuck2.snr = zeros(length(subsampling),numTrials);
tuck2.times = zeros(length(subsampling),numTrials);
tuck3 = struct;
tuck3.snr = zeros(length(subsampling),numTrials);
tuck3.times = zeros(length(subsampling),numTrials);
htuck1 = struct;
htuck1.snr = zeros(length(subsampling),numTrials);
htuck1.times = zeros(length(subsampling),numTrials);
htuck2 = struct;
htuck2.snr = zeros(length(subsampling),numTrials);
htuck2.times = zeros(length(subsampling),numTrials);
htuck3 = struct;
htuck3.snr = zeros(length(subsampling),numTrials);
htuck3.times = zeros(length(subsampling),numTrials);

maxIter = 100;
verbosity = 0;

%Save file + log file name
saveName = 'singlereflector-missingall';

%Log file
fid = fopen([saveFigsDir saveName '.txt'],'w');
if exist('metadata','var')
    fprintf(fid,['Metadata: ' metadata '\n']);
end
fprintf(fid,['-----------------------------------------------\n']);
fprintf(fid,['T1 - tucker interpolation with ranks ' num2str(tuckRanks1) '\n']);
fprintf(fid,['T2 - tucker interpolation with ranks ' num2str(tuckRanks2) '\n']);
fprintf(fid,['T3 - tucker interpolation with ranks ' num2str(tuckRanks2) ', rank continuation\n']);
fprintf(fid,['HT1 - hierarchical tucker interpolation with ranks kleaf = ' num2str(kleaf1) ' , kint = ' num2str(kint1) ' \n']);
fprintf(fid,['HT2 - hierarchical tucker interpolation with ranks kleaf = ' num2str(kleaf2) ', kint = ' num2str(kint2) ' \n']);
fprintf(fid,['HT3 - hierarchical tucker interpolation with ranks kleaf = ' num2str(kleaf3) ', kint = ' num2str(kint3) ' \n']);
fprintf(fid,['-----------------------------------------------\n']);


for i=1:length(subsampling)    
    disp(['Experiment ' num2str(i)]);
    for j=1:numTrials
        disp(['Trial ' num2str(j) ' of ' num2str(numTrials)]);
        % Generate test, train data               
        trainSet = randperm(prod(dims),round(subsampling(i)*prod(dims)));
        Rtrain = opRestriction(prod(dims),trainSet);
        Rtest = opRestriction(prod(dims),setdiff(1:prod(dims),trainSet));
        
        trainSet = logical( Rtrain' * ones(size(Rtrain,1),1));        
        trainData = fslice;
        trainData(~trainSet) = 0;
        testData = Rtest * vec(fslice);                

        % Tuck1 experiment
        expTime = tic;
        [U,B] = fitTucker(trainSet, trainData, dims,tuckRanks1,'verbosity',verbosity,'maxIter',maxIter); 
        expTime = toc(expTime);
                        
        X_tuck1 = ttm(B,U,1:4);
        
        snr_t1 = SNR(testData,Rtest*vec(ttm(B,U,1:4)));
        tuck1.snr(i,j) = snr_t1;
        tuck1.times(i,j) = expTime;
        
        % Tuck2 experiment
        expTime = tic;
        [U,B] = fitTucker(trainSet, trainData, dims,tuckRanks2,'verbosity',verbosity,'maxIter',maxIter); 
        expTime = toc(expTime);
                
        X_tuck2 = ttm(B,U,1:4);
        snr_t2 = SNR(testData,Rtest*vec(ttm(B,U,1:4)));
        tuck2.snr(i,j) = snr_t2;
        tuck2.times(i,j) = expTime;
        
        % Tuck3 experiment       
        expTime = tic;
        [U,B] = fitTuckerRankIncrease(trainSet, trainData, dims, tuckRanks2, 5, 'verbosity',verbosity,'maxIter',maxIter);
        expTime = toc(expTime);
                
        snr_t3 = SNR(testData,Rtest*vec(ttm(B,U,1:4)));
        tuck3.snr(i,j) = snr_t3;   
        tuck3.times(i,j) = expTime;
        X_tuck3 = ttm(B,U,1:4);
        
        %HTuck1 experiment
        expTime = tic;
        x = fitHT(trainSet,trainData,dimTree1,'verbosity',verbosity,'progTol',1e-6,'maxIter',maxIter);
        expTime = toc(expTime);
        
        X_htuck1 = dimTree1.fullND(x);
        
        snr_h1 = SNR(testData,Rtest*(dimTree1.full(x)));
        htuck1.snr(i,j) = snr_h1;
        htuck1.times(i,j) = expTime;
        
        %HTuck2 experiment
        expTime = tic;
        x = fitHT(trainSet,trainData,dimTree2,'verbosity',verbosity,'progTol',1e-6,'maxIter',maxIter);
        expTime = toc(expTime);
        
        X_htuck2 = dimTree2.fullND(x);
        
        snr_h2 = SNR(testData,Rtest*(dimTree2.full(x)));
        htuck2.snr(i,j) = snr_h2;
        htuck2.times(i,j) = expTime;
        
        %HTuck3 experiment
        expTime = tic;
        x = fitHT(trainSet,trainData,dimTree3,'verbosity',verbosity,'progTol',1e-6,'maxIter',maxIter);
        expTime = toc(expTime);
        
        snr_h3 = SNR(testData,Rtest*(dimTree3.full(x)));
        htuck3.snr(i,j) = snr_h3;
        htuck3.times(i,j) = expTime; 
        
        X_htuck3 = dimTree3.fullND(x);
        
        save([saveDir 'singlereflector-missingall.mat'],'tuck1','tuck2','tuck3','htuck1','htuck2','htuck3');
        %Output a common source gather for a particular subsampling
        %ratio
         if i==1 && j==1 || i==2 && j==1
            visible = true;
            srcx = 20; srcy = 20;
            slice = @(A) squeeze(A(:,srcx,:,srcy));
            true_slice = slice(fslice);
            imagePlot(true_slice,'cbar',true,'visible',visible,'centercaxis',true);

            sub_str = ['sub-' num2str(subsampling(i)) '-'];
            saveFigure = @(filename,postfilename) print(gcf,'-depsc',[saveFigsDir sub_str filename '-' postfilename '.eps']);
            
            saveFigure('true','');
            
            true_axis = caxis;
            
            interpPlot = @(A) imagePlot(slice(A),'cbar',true,'visible',visible,'coloraxis',true_axis);
                        
            interpPlot(trainData);
            saveFigure('trainData','');
            
            
            interpPlot(X_tuck2);           
            saveFigure('interp-tuck2',['snr-' num2str(SNR(vec(true_slice),vec(slice(X_tuck2))),'%3.3e')]);
            
            interpPlot(X_tuck1);            
            saveFigure('interp-tuck1',['snr-' num2str(SNR(vec(true_slice),vec(slice(X_tuck1))),'%3.3e')]);
            
            interpPlot(X_tuck3);
            saveFigure('interp-tuck3',['snr-' num2str(SNR(vec(true_slice),vec(slice(X_tuck3))),'%3.3e')]);
            
            interpPlot(X_htuck3);            
            saveFigure('interp-htuck3',['snr-' num2str(SNR(vec(true_slice),vec(slice(X_htuck3))),'%3.3e')]);
            
            interpPlot(X_htuck2);
            saveFigure('interp-htuck2',['snr-' num2str(SNR(vec(true_slice),vec(slice(X_htuck2))),'%3.3e')]);   
            
            interpPlot(X_htuck1);
            saveFigure('interp-htuck1',['snr-' num2str(SNR(vec(true_slice),vec(slice(X_htuck1))),'%3.3e')]);  
            
        end        
        
    end    

end
% Save experiment results
save([saveDir 'singlereflector-missingall.mat'],'tuck1','tuck2','tuck3','htuck1','htuck2','htuck3');

%% Experiment output
fprintf(fid,'Recovery\n');
disp('Recovery');
f = @(x) mean(x,2);
mean_snr = [f(tuck1.snr), f(tuck2.snr), f(tuck3.snr), f(htuck1.snr), f(htuck2.snr), f(htuck3.snr)];
f = @(x) std(x,0,2);
std_snr = [f(tuck1.snr), f(tuck2.snr), f(tuck3.snr), f(htuck1.snr), f(htuck2.snr), f(htuck3.snr)];

sigfigs = 2;
results = matrix2table(mean_snr,std_snr,{'Percent receivers kept', 'T1','T2','T3','HT1','HT2','HT3'},arrayfun(@num2str,subsampling,'UniformOutput',false),sigfigs);
for i=1:length(results)        
    fprintf(fid,[results{i} '\n']);
    disp(results{i});
    if i < length(results)
        fprintf(fid,[repmat('-',1,length(results{i})) '\n']);
        disp(repmat('-',1,length(results{i})));
    end
end


fprintf(fid,'Times\n');
disp('Times');
f = @(x) mean(x,2);
mean_times = [f(tuck1.times), f(tuck2.times), f(tuck3.times), f(htuck1.times), f(htuck2.times),f(htuck3.times)];
f = @(x) std(x,0,2);
std_times = [f(tuck1.times), f(tuck2.times), f(tuck3.times), f(htuck1.times), f(htuck2.times),f(htuck3.times)];
sigfigs = 2;
results = matrix2table(mean_times,std_times,{'Percent data kept', 'T1','T2','T3','HT1','HT2','HT3'},arrayfun(@num2str,subsampling,'UniformOutput',false),sigfigs);
for i=1:length(results)
    fprintf(fid,[results{i} '\n']);
    disp(results{i});
    if i < length(results)
        fprintf(fid,[repmat('-',1,length(results{i})) '\n']);
        disp(repmat('-',1,length(results{i})));
    end
end

fclose(fid);

%% Output results to latex, for the paper
fid = fopen([saveFigsDir 'latex.txt'],'w');
snr_format = '%3.1f';
time_format = '%5d';
fprintf(fid, ['\\begin{table}\n \\centering \n' ...
              '\\begin{tabular}{c | c c | c c | c c } \n' ...
              '\\Xhline{2\\arrayrulewidth}  \\\\ \n' ...
              '& \\multicolumn{6}{c}{\\centering Single reflector ' ...
              'data - sampling percentage (missing points) }\\\\ \n' ...
              '&\\multicolumn{2}{c}{\\centering 10\\%%} & \\multicolumn{2}{c}{\\centering 30\\%%} & \\multicolumn{2}{c}{\\centering 50\\%%} \\\\ \n' ...
              '& SNR [dB] & time [s] & SNR [dB] & time [s] & SNR [dB] & time [s] \\\\ \n' ...
              '\\hline\n'] );

fprintf(fid, ['geomCG(' num2str(tuckRanks1(1)) ') - 0 & ' num2str(mean_snr(1,1),snr_format) ' & ' num2str(round(mean_times(1,1)),time_format) ' & ' num2str(mean_snr(2,1),snr_format) ' & ' num2str(round(mean_times(2,1)),time_format) ' & ' num2str(mean_snr(3,1),snr_format) ' & ' num2str(round(mean_times(3,1)),time_format) ' \\\\ \n']);

fprintf(fid, ['geomCG(' num2str(tuckRanks2(1)) ') - 0 & ' num2str(mean_snr(1,2),snr_format) ' & ' num2str(round(mean_times(1,2)),time_format) ' & ' num2str(mean_snr(2,2),snr_format) ' & ' num2str(round(mean_times(2,2)),time_format) ' & ' num2str(mean_snr(3,2),snr_format) ' & ' num2str(round(mean_times(3,2)),time_format) ' \\\\ \n']);

fprintf(fid, ['geomCG(' num2str(tuckRanks2(1)) ') - 5 & ' num2str(mean_snr(1,3),snr_format) ' & ' num2str(round(mean_times(1,3)),time_format) ' & ' num2str(mean_snr(2,3),snr_format) ' & ' num2str(round(mean_times(2,3)),time_format) ' & ' num2str(mean_snr(3,3),snr_format) ' & ' num2str(round(mean_times(3,3)),time_format) ' \\\\ \n']);

fprintf(fid, ['HTOpt(' num2str(kleaf1) ',' num2str(kint1) ') & ' num2str(mean_snr(1,4),snr_format) ' & ' num2str(round(mean_times(1,4)),time_format) ' & ' num2str(mean_snr(2,4),snr_format) ' & ' num2str(round(mean_times(2,4)),time_format) ' & ' num2str(mean_snr(3,4),snr_format) ' & ' num2str(round(mean_times(3,4)),time_format) ' \\\\ \n']);

fprintf(fid, ['HTOpt(' num2str(kleaf2) ',' num2str(kint2) ') & ' num2str(mean_snr(1,5),snr_format) ' & ' num2str(round(mean_times(1,5)),time_format) ' & ' num2str(mean_snr(2,5),snr_format) ' & ' num2str(round(mean_times(2,5)),time_format) ' & ' num2str(mean_snr(3,5),snr_format) ' & ' num2str(round(mean_times(3,5)),time_format) ' \\\\ \n']);

fprintf(fid, ['HTOpt(' num2str(kleaf3) ',' num2str(kint3) ') & ' num2str(mean_snr(1,6),snr_format) ' & ' num2str(round(mean_times(1,6)),time_format) ' & ' num2str(mean_snr(2,6),snr_format) ' & ' num2str(round(mean_times(2,6)),time_format) ' & ' num2str(mean_snr(3,6),snr_format) ' & ' num2str(round(mean_times(3,6)),time_format) ' \\\\ \n']);


fprintf(fid,['\\Xhline{2\\arrayrulewidth}\n'...
             '\\end{tabular}\n'...
             '\\caption{Reconstruction results for single reflector data - missing points - mean SNR over ' num2str(numTrials) ' random training test sets} \n'...
             '\\label{table:singlereflector-missingall} \n' ...
             '\\end{table} \n']);
fclose(fid);







