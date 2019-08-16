function [x,dimTree,lambda] = crossValidate(rhs, e, dimTree, ...
                                            minRanks,...
                                            maxRanks,...
                                            minLambda,...
                                            maxLambda)

if exist('minLambda','var')== 0
    minLambda = [];
end
if exist('maxLambda','var')==0
    maxLambda = [];
end

numSteps = 5;
nIter = 25;
method = 'GN';

trainPercent = 0.8;
trainI = randperm(length(find(e)),round(trainPercent*length(find(e))));
trainSet = e;
trainSet(~trainI) = 0;
testSet = e & ~trainSet;

dRanks = cell_skeleton(minRanks);

nRankIncreases = 0;
for i=1:length(dRanks)
    for j=1:length(dRanks{i})
        dRanks{i}{j} = floor((maxRanks{i}{j} - minRanks{i}{j})/numSteps);
        nRankIncreases = nRankIncreases + 1;
    end
end
    
if ~isempty(maxLambda) && minLambda > 0 && ~isempty(maxLambda) && maxLambda > 0
    dLambda = [0,logspace(minLambda,maxLambda,numSteps-1)];
else
    dLambda = [0];
end


dimTreeOld = dimTree;

for i=1:nRankIncreases*numSteps
    
end





end