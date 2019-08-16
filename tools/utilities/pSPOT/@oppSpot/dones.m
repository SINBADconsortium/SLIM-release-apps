function y = dones(op,varargin)
%DONES Distributed ones vector in operator domain
%
%   y = dones(A) generates a ones vector with the size of the operator
%   domain, distributed accordingly to the operators needs so that A*y is
%   a valid operation.
%
%   y = dones(A,NCOLS) generates NCOLS vectors

if isempty(varargin)
    ncols = 1;
else
    ncols = varargin{1};
end

if length(op.opsn) > 1 % Distributed
    scheme  = pSPOT.utils.compositeDef(op.opsn);
    glosize = [op.n ncols];

    spmd
        scheme  = sum(scheme);
        ypart   = ones(scheme,ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme;
        ycodist = codistributor1d(1,ygpart,glosize);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
    
else % Non-distributed
    y = ones(op.n,ncols);
end