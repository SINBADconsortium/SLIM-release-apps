function y = rrandn(op,varargin)
%RRANDN Distributed random vector in operator range
%
%   y = rrandn(A) generates a random vector y with the size of the operator
%   range, distributed accordingly to the operators needs so that A'*y is
%   a valid operation.
%
%   y = rrandn(A,NCOLS) generates NCOLS vectors

ncols = [varargin{:}];
if isempty(ncols)
    ncols = 1;
end

if length(op.opsm) > 1 % Distributed
    scheme  = pSPOT.utils.compositeDef(op.opsm);
    glosize = [op.m ncols];

    spmd
        scheme  = sum(scheme);
        ypart   = randn(scheme,ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme;
        ycodist = codistributor1d(1,ygpart,glosize);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
    
else % Non-distributed
    y = randn(op.m,ncols);
end