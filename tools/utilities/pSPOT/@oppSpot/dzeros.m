function y = dzeros(op,varargin)
%DZEROS Distributed zero vector in operator domain
%
%   y = dzeros(A) generates a zero vector with the size of the operator
%   domain, distributed accordingly to the operators needs so that A*y is
%   a valid operation.
%
%   y = drandn(A,NCOLS) generates NCOLS vectors

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
        ypart   = zeros(scheme,ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme;
        ycodist = codistributor1d(1,ygpart,glosize);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
    
else % Non-distributed
    y = zeros(op.n,ncols);
end