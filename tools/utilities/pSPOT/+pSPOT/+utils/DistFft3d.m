function X = DistFft3d( X )
%DISTFFT3D(X) takes in a distributed nd-array of the arrangement X = [T X Y]
%and runs fft along each dim, then returns an nd-array arranged as X =
%[ky kx f]

    spmd
        codist = getCodistributor(X);
        distdim = codist.Dimension;
        part = codist.Partition;
        sizes = size(X);
        
        %fft along T
        if 1 == distdim
            X = fft(X);
            X = getLocalPart(X);
        else
            X = getLocalPart(X);
            X = fft(X);
        end
        
        %fft along X
        if 2 == distdim
            X = codistributed.build(X, codistributor1d(distdim,part,sizes));
            X = redistribute(X,codistributor1d(1));
            X = getLocalPart(X);
            part = codistributor1d.unsetPartition;
        end
        X = permute(X, [2 1 3]);
        X = fft(X);
        sizes = [sizes(2), sizes(1), sizes(3)];
        
        %fft along Y
        if 3 == distdim
            X = codistributed.build(X, codistributor1d(distdim,part,sizes));
            X = redistribute(X,codistributor1d(2));
            X = getLocalPart(X);
            part = codistributor1d.unsetPartition;
        end
        X = permute(X, [3 1 2]);
        X = fft(X);
        sizes = [sizes(3), sizes(1), sizes(2)];

        %return volume as [ky kx f] distributed along f
        X = codistributed.build(X,codistributor1d(3,part,sizes));
    end
    
end             