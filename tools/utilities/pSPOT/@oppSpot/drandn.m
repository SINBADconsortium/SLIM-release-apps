function y = drandn(op,varargin)
%DRANDN Distributed random vector in operator domain
%
%   y = drandn(A) generates a random vector y with the size of the operator
%   domain, distributed accordingly to the operators needs so that A*y is
%   a valid operation.
%
%   y = drandn(A,NCOLS) generates NCOLS vectors

ncols = [varargin{:}];
if isempty(ncols)
    ncols = 1;
end

if length(op.opsn) > 1 % Distributed
    scheme  = pSPOT.utils.compositeDef(op.opsn);
    glosize = [op.n ncols];

    spmd
        scheme  = sum(scheme);
        ypart   = randn(scheme,ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme;
        ycodist = codistributor1d(1,ygpart,glosize);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
    
else % Non-distributed
    y = randn(op.n,ncols);
end

