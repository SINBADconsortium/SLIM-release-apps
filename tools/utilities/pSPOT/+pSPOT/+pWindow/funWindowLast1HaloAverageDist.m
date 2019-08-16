function Y = funWindowLast1HaloAverageDist( x, n, p, h )
%funWindowLast1HaloAverageDist
%       is a support function for forward opdWindowLast1* operators;
%       it exchanges halos of the input vector X and avarages ovelaps
%
%   Y = funWindowLast1HaloAverageDist( X, N, P, H )
%
%   INPUT:
%      X = input vector
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      Y = output vector

    assert(isvector(x),'Fatal error: x has to be vector')
    assert(isdistributed(x),'Fatal error: x has to be distributed')
    assert(parpool_size()==p,'Fatal error: p does not match parallel pool size')
    assert(h>0,'Fatal error: half-halo has to be beger than 0')
    N=n+(p-1)*2*h;
    mN=prod(size(x));
    assert(mod(mN,N)==0,'Fatal error: N is not a valid last dimension')
    m=prod(size(x))/N;
    spmd
        part=codistributor1d.defaultPartition(n);
        PART=zeros(1,p);
        PART(1)=part(1)+h;
        PART(2:end-1)=part(2:end-1)+2*h;
        PART(end)=part(end)+h;
    
        % get local part
        mydata = getLocalPart(x);
        % test m*PART(labindex) here
        assert(prod(size(mydata))==m*PART(labindex),'Fatal error: wrong local part: partition is not default?')
        mydata = reshape(mydata,[m PART(labindex)]);
        otherdata=mydata;
    
        % shift right
        if labindex<numlabs; labTo = labindex + 1; else labTo = []; end;
        if labindex>1; labFrom = labindex - 1; else labFrom = []; end;
        %fprintf('shift right: %d <%d >%d\n',labindex,labFrom,labTo);
        halo = labSendReceive(labTo, labFrom, mydata(:,end-2*h+1:end));
        if labindex > 1; otherdata(:,1:1+2*h-1)=(otherdata(:,1:1+2*h-1)+halo)/2.; end;
    
        % shift left
        if labindex>1; labTo = labindex - 1; else labTo = []; end;
        if labindex<numlabs; labFrom = labindex + 1; else labFrom = []; end;
        %fprintf('shift left: %d <%d >%d\n',labindex,labFrom,labTo);
        halo = labSendReceive(labTo, labFrom, mydata(:,1:1+2*h-1));
        if labindex < numlabs; otherdata(:,end-2*h+1:end)=(otherdata(:,end-2*h+1:end)+halo)/2.; end;
    
        % rebuild array
        codist = codistributor1d(1,m*PART,[m*N 1]);
        Y=codistributed.build(otherdata(:),codist,'noCommunication');
    
    end
end
