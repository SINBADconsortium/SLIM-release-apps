function Y = funWindowLast1HaloDropCodist( x, n, p, h )
%funWindowLast1HaloDropCodist
%       is a support function for forward opdWindowLast1* operators;
%       it drops halos from the input vector X
%
%   Y = funWindowLast1HaloDropCodist( X, N, P, H )
%
%   INPUT:
%      X = input vector
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      Y = output vector

    assert(isvector(x),'Fatal error: x has to be vector')
    assert(iscodistributed(x),'Fatal error: x has to be codistributed')
    assert(numlabs==p,'Fatal error: p does not match parallel pool size')
    assert(h>0,'Fatal error: half-halo has to be beger than 0')
    N=n+(p-1)*2*h;
    mN=prod(size(x));
    assert(mod(mN,N)==0,'Fatal error: N is not a valid last dimension')
    m=prod(size(x))/N;
    %spmd
        part=codistributor1d.defaultPartition(n);
        PART=zeros(1,p);
        PART(1)=part(1)+h;
        PART(2:end-1)=part(2:end-1)+2*h;
        PART(end)=part(end)+h;
    
        % get local part
        mydata = getLocalPart(x);
        % test m*PART(labindex) here
        assert(prod(size(mydata))==m*PART(labindex),'Fatal error: wrong local part: partition is not the default?')
        mydata = reshape(mydata,[m PART(labindex)]);
        if labindex~=1 &  labindex~=numlabs
            otherdata=mydata(:,1+h:end-h);
        elseif labindex==1
            otherdata=mydata(:,1:end-h);
        elseif labindex==numlabs
            otherdata=mydata(:,1+h:end);
        end
    
        % rebuild array
        codist = codistributor1d(1,m*part,[m*n 1]);
        Y=codistributed.build(otherdata(:),codist,'noCommunication');
    
    %end
end
