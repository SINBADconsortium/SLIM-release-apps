function [U,B] = fitTuckerRankIncrease( trainIndices, rhs, dims, ...
                                        ranks, numRankIncreases, varargin )
% FITTUCKERRANKINCREASE - Incremental rank increase strategy for interpolation of Tucker tensors,
% using the same strategy as 'Low-Rank Tensor Completion by Riemannian Optimization'
% by D. Kressner, M. Steinlechner, B Vandereycken.
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca 
%
% Usage:
%    [U,B] = fitTuckerRankIncrease( trainIndices, rhs, dims, ranks,
%                                   numRankIncreases, ... )
% Input: 
%    trainIndices     - 1/0 logical vector of length == length(rhs), which should satisfy
%                       trainIndices(i) == true if and only if rhs(i) is known
%    rhs              - zero-filled training data
%    dims             - dimensions of the target tensor
%    ranks            - Tucker ranks of the target tensor (need length(ranks) == length(dims) )
%    numRankIncreases - number of times to increase ranks from
%                       ones(1,length(U)) to the target rank
% 
% Optional Input:
%    All of the same optional inputs as fitTucker.m, refer to documentation
%    there for more details
                                    
d = length(dims);
function v = vecParams(U,B)
v = [];
for j=1:length(U)
    v = [v;vec(U{j})];
end
v = [v;vec(B)];
end

r = ones(1,d);
U = cell(d,1);
for i=1:d
    U{i} = orth(randn(dims(i),1));
end
B = randn(r);    

x0 = vecParams(U,B);

dRanks = floor((ranks-1)/numRankIncreases);
dRanks = max(dRanks,1);

for i=0:numRankIncreases-1
    if i==numRankIncreases-1
        args = varargin;
        args{end+1} = 'progTol'; args{end+1} = 1e0;
        args{end+1} = 'x0'; args{end+1} = x0;
        [U,B] = fitTucker(trainIndices, rhs, dims, ranks, args{:});
    else
        curRanks = max(min(r + i * dRanks,ranks),1);
        args = varargin;
        args{end+1} = 'progTol'; args{end+1} = 1e-3;
        args{end+1} = 'x0'; args{end+1} = x0;
        [U,B] = fitTucker(trainIndices, rhs, dims,curRanks ,args{:});    
        
        if i==numRankIncreases-2
            newRanks = ranks;
        else
            newRanks = min(r + (i+1)*dRanks,ranks);
        end
        Bnew = zeros(newRanks);
        subIdx = '';
        for i=1:d
            subIdx = [subIdx,'1:size(B,' num2str(i) '),'];
        end
        subIdx = subIdx(1:end-1);
        eval(['Bnew(' subIdx ') = B;']);

        B = Bnew;
        for i=1:d
            Ui = zeros(dims(i),newRanks(i));
            Ui(:,1:curRanks(i)) = U{i};
            U{i} = Ui;
        end
        x0 = vecParams(U,B);
    end        
end


end