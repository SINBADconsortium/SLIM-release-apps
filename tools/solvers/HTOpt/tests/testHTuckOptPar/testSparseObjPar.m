function testSparseObjPar
% Test that the dense + sparse objectives are giving the same result
% In parallel
    n = 10; d = 5; dims = n*ones(1,d);
    kint = 10; kleaf = 5;
    dimTree = dimensionTree(dims,kleaf,kint);
    x = project(dimTree.randn(),dimTree);
    numSamplePts = 100*n;
    
    e = distributed.false(prod(dims),1);
    
    spmd
        eloc = getLocalPart(e);
        cdist = getCodistributor(e);
        if mod(numSamplePts,numlabs)==0
            samplePtsPerLab = numSamplePts/numlabs;
        else
            samplePtsPerLab = floor(numSamplePts/numlabs);
            if labindex == numlabs
                samplePtsPerLab = samplePtsPerLab + ...
                    mod(numSamplePts,numlabs);
            end
        end        
        
        idx = globalIndices(e,1);
        minIdx = min(idx); maxIdx = max(idx);
        
        I = minIdx + sort(randperm(maxIdx-minIdx+1,samplePtsPerLab))-1;
        
        eloc(I+1-minIdx) = true;
        
        e = codistributed.build(eloc,cdist);
        Ifullloc = cell(1,d);
        [Ifullloc{:}] = ind2sub(dims,I);
        Ifullloc = reshape(cell2mat(Ifullloc),length(I),d);

        partSize = floor(numSamplePts/numlabs) * ones(1,numlabs);
        partSize(end) = partSize(end) + mod(numSamplePts,numlabs);
        
        codist = codistributor1d(1,partSize,[numSamplePts,d]);
        Ifull = codistributed.build(Ifullloc,codist,'noCommunication'); 
    end    
    spmd
        partSize = getCodistributor(Ifull);
        partSize = partSize.Partition;
        codist = codistributor1d(1,partSize,[size(Ifull,1),1]);
        bsub = codistributed.randn(size(Ifull,1),1,codist);    
        codist = getCodistributor(e);   
        bsub_full = codistributed.zeros(size(e,1),1,codist);
    end    
    bsub_full(find(e)) = bsub;
    Ifull = Ifull';
    %Dense, parallel code

    X = dimTree.fullDist(x);
    rfull = pSPOT.utils.distVectorize(X);
    spmd
        rfullloc = getLocalPart(rfull);
        eloc = getLocalPart(e);
        rfullloc(~eloc) = 0;
        codist = getCodistributor(rfull);
        rfull = codistributed.build(rfullloc,codist,'noCommunication');
    end
    rfull = (rfull - bsub_full);
    f = 0.5*norm(rfull,2)^2;
    J = oppHTuckJ2(dimTree,x);
    g = project_horizontal(x,J' * rfull,dimTree);
    
    %Sparse, parallel code
    [f2,g2] = LSMisfitHT(x,Ifull,bsub,dimTree);
    
    assertVectorsAlmostEqual(f,f2);
    assertVectorsAlmostEqual(g,g2);
           

