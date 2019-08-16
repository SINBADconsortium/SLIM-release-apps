function [fval,dL,dR,r] = spgLRobj(L,R,b,e,penalty)
% SPGLROBJ - SPGLR objective function
%    Computes    
%       sum_{i,j} penalty( L(i,:)'*R(j,:) - b(i,j) )
%    where the sum runs over all i,j, such that e(i,j) is true and
%    L, R,b are given by distributed arrays. Penalty is any
%    differentiable function. 
%
%
% Usage:
%   [fval,dL,dR] = spgLRobj(L,R,b,e,penalty)
% 
% Input:
%   L,R     - L,R blocks distributed over rows
%   b       - incomplete data matrix, distributed identically to L,R
%   e       - 0-1 matrix, same size and distribution as b, 1 where
%             there is a data point, 0 otherwise
%   penalty - differentiable function handle, applied pointwise to
%             residual matrix
    
    gradMode = nargout > 1; computeRes = nargout >= 4; 
    nlabs = parpool_size();
    nrows = size(b,1); ncols = size(b,2);
    assert( nlabs > 0, 'Need open parallel pool');
    rowidx = ceil(linspace(1,nrows+1,nlabs+1));
    colidx = ceil(linspace(1,ncols+1,nlabs+1));
    rank = size(L,2);
    spmd                        
        Lloc = getLocalPart(L); Rloc = getLocalPart(R);          
        eloc = getLocalPart(e); dataLoc = getLocalPart(b);
        floc = 0;  
        if gradMode
            dLloc = zeros(size(Lloc)); dRloc = zeros(size(Rloc));
            recvblk = [vec(Rloc); vec(dRloc)];
        else
            dLloc = []; dRloc = [];
            recvblk = Rloc; 
        end
        T = circshift(vec(1:nlabs),-(labindex-1))';
        if computeRes
            ixloc = zeros(length(find(eloc)),1); 
            iyloc = zeros(length(find(eloc)),1);
            rblkloc = zeros(length(find(eloc)),1);                
            rel = 1;
        end
        sendIdx = mod(labindex,nlabs)+1; recvIdx = mod(labindex-2,nlabs)+1;      
        % iterate over cyclic permutation of the blocks
        for ic=1:length(T)
            nextIdx = mod(T(ic)+1,nlabs+1); if nextIdx==0,nextIdx=1;end
            minCol = colidx(T(ic)); if nextIdx==1,maxCol = ncols; else maxCol = colidx(nextIdx)-1;end                               
            if gradMode
                Rloc = reshape(recvblk(1:length(recvblk)/2),maxCol-minCol+1,rank);
                dRloc = reshape(recvblk(length(recvblk)/2+1:end),maxCol-minCol+1,rank);
            else
                Rloc = recvblk;
            end            
            % Get data corresponding to this block                
            bloc = dataLoc(:,minCol:maxCol);
            elocblk = eloc(:,minCol:maxCol);
            
            % Local computation
            rloc = Lloc * Rloc';
            rloc(~elocblk) = 0;
            rloc(elocblk) = rloc(elocblk)-bloc(elocblk);                
            [floc1,rloc] = penalty(rloc);
            floc = floc + floc1;
            
            if computeRes                                 
                [ixblk,iyblk] = ind2sub([size(Lloc,1),size(Rloc,1)],find(elocblk));                   
                ixloc(rel:length(ixblk)+rel-1) = ixblk+minRow-1;
                iyloc(rel:length(iyblk)+rel-1) = iyblk+minCol-1;
                rloc = rloc(elocblk);
                rblkloc(rel:rel+length(rloc)-1) = rloc;
                rel = rel+length(iyblk);
            end

            % Communicate R block, dR block to next worker if necessary
            if gradMode
                dLloc = dLloc + rloc * Rloc;
                dRloc = dRloc + rloc' * Lloc;
                sendBlk = [vec(Rloc); vec(dRloc)];
            else
                sendBlk = Rloc; 
            end
            
            recvblk = labSendReceive(recvIdx,sendIdx,sendBlk);
        end
        if gradMode
            Rloc = reshape(recvblk(1:length(recvblk)/2),length(recvblk)/(2*rank),rank);
            dRloc = reshape(recvblk(length(recvblk)/2+1:end),length(recvblk)/(2*rank),rank);            
            dL = codistributed.build(dLloc,codistributor1d(1,diff(rowidx),[nrows,rank]),'noCommunication');
            dR = codistributed.build(dRloc,codistributor1d(1,diff(colidx),[ncols,rank]),'noCommunication');
        end
        if computeRes   
            ixloc = codistributed.build(ixloc,codistributor1d(1,nnz',[sum(nnz),1]),'noCommunication');
            iyloc = codistributed.build(iyloc,getCodistributor(ixloc),'noCommunication');
            rblk = codistributed.build(rblkloc,getCodistributor(ixloc),'noCommunication');
            r = sparse(ixloc,iyloc,rblk,nrows,ncols);
        end

        fval = gplus(floc,1);
    end
    fval = sqrt(fval{1});
    if gradMode
        dL = dL*(1/(2*fval)); dR = dR*(1/(2*fval)); 
    end

end
