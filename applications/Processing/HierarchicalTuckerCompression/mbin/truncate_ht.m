function [dimTree,x] = truncate_ht(X, err, relerr,randsvd)
% TRUNCATE_HT - Truncates a full tensor to HT format, given a user-supplied error
% bound. Returns a new dimension tree + set of HT paramters.
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
%
% Use:
%  [dimTree,x] = truncate_ht(X,{err},{relerr},{randsvd})
%
% Input:
%  X       - full tensor
%  err     - error tolerance (default: 1e-3)
%  relerr  - true (default) - relative error bound
%          - false - absolute error bound
%  randsvd - if true, uses a randomized SVD instead of a regular SVD
%            on the full matricization (default: false)
%
% Output:
%  dimTree - dimension tree object with the properly selected ranks
%  x       - vectorized HT parameters

if exist('err','var')==0
    err = 1e-3;
end
if exist('relerr','var')==0
    relerr = true;
end
if exist('randsvd','var')==0
    randsvd = false;
end

    function k = optimal_rank(svals)
        delta = zeros(length(svals));
        for i=1:length(delta)
            delta(i) = sum(svals(i+1:end).^2);
        end
        k = find(sqrt(delta) <= err_bound, 1 );
    end


dimTree = dimensionTree(size(X),1,1);
[U,B] = dimTree.fromVec(dimTree.randn());
T = dimTree.emptyTree(U,B);

if relerr
    err_bound = norm(vec(X)) * err/sqrt(2*length(size(X))-3);
else
    err_bound = err/sqrt(2*length(size(X))-3);
end


leafRank = 0;
itr = dimTree.iterator('leaves',T);
C_lower = X;
while itr.advance()
    dim = itr.getDims();
    xmat = matricize(X,dim);
    [leftSV,S] = svd(xmat * xmat');
    k = optimal_rank(sqrt(diag(S)));
    U = leftSV(:,1:k); itr.setValue('U',U);
    itr.setValue('k',k);
    leafRank = max(leafRank,k);
    C_lower = ttm(C_lower, U',dim);
    clear leftSV S U;
    
end
C_level = C_lower;
dimPartition = {};
d = 1;
for i=1:length(B{end})
    height = length(B);
    if dimTree.isLeaf(height,i)
        dimPartition{end+1} = [d];
        d = d + 1;
    else
        dimPartition{end+1} = [d,d+1];
        d = d + 2;
    end
end
dimOffset = 1;
itr = dimTree.iterator('up',itr.T);
while itr.advance()
    if ~itr.isRoot()
        if ~itr.isLeaf()
            if randsvd
                [leftSV,S,~] = randsvd(matricize(C_level,dimPartition{dimOffset}),500,1e-5);
            else
                [leftSV,S,~] = svd(matricize(C_level,dimPartition{dimOffset}),0);
            end
            k = optimal_rank(diag(S));
            B = dematricize(leftSV(:,1:k),[itr.getLeftValue('k'),itr.getRightValue('k'),k],[1 2]);
            itr.setValue('B',B);
            itr.setValue('k',k);
            clear leftSV S;
            
            %Tensor contraction
            beforeDims = 1:dimOffset-1;
            afterDims = dimOffset+2:length(size(C_lower));
            
            C_lower = ttt(B,C_lower,[1 2],dimOffset:dimOffset+1);
            % Consistent dimension ordering
            if ~isempty(afterDims)
                afterDims = afterDims-1;
            end
            if ~isempty(beforeDims)
                beforeDims = beforeDims+1;
            end
            C_lower = permute(C_lower,[beforeDims 1 afterDims]);
        end
        dimOffset = dimOffset + 1;
        if itr.idx == length(dimTree.nodes{itr.depth})
            level = itr.depth;
            C_level = C_lower;
            dimOffset = 1;
            if level > 2
                dimPartition = {};
                d = 1;
                for i=1:length(dimTree.nodes{level-1})
                    if dimTree.isLeaf(level-1,i)
                        dimPartition{end+1} = [d];
                        d = d + 1;
                    else
                        dimPartition{end+1} = [d,d+1];
                        d = d + 2;
                    end
                end
            end
        end
    else
        B = C_level;
        itr.setValue('B',B);
    end
end
[U,B] = dimTree.extractParams(itr.T,'U','B');

itr = dimTree.iterator('down',itr.T);
while itr.advance();
    if ~itr.isRoot()
        dimTree.rank(itr.depth,itr.idx,itr.getValue('k'));
    end
end
x = dimTree.toVec(U,B);
end