function X = DistFreqSlice(X)
%DISTFREQSLICE(X) takes in a T,X,Y,.. distributed nd-array and returns a 
%  X,Y,..,f freq slice volume distributed by the frequencies.

assert( isa(X,'distributed'), 'input must be distributed');

   sizes = size(X);
   perm = [ 2:ndims(X), 1];
   
   spmd
     codist = getCodistributor(X);
     part = codist.Partition;
     dim = codist.Dimension;
     X = fft(X);
     X = getLocalPart(X);
     X = permute(X,perm);
     sizes = circshift(sizes,[0 -1]);
     X = codistributed.build(X,codistributor1d(dim-1,part,sizes));
     X = redistribute(X,codistributor1d(ndims(X)));
   end

end