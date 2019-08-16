function y = rones(op,varargin)
%RONES Distributed ones vector in operator range
%
%   y = rones(A) generates a ones vector with the size of the operator
%   range, distributed accordingly to the operators needs so that A'*y is
%   a valid operation.
%
%   y = rones(A,NCOLS) generates NCOLS vectors

if isempty(varargin)
    ncols = 1;
else
    ncols = varargin{1};
end

if length(op.opsm) > 1 % Distributed
    scheme  = pSPOT.utils.compositeDef(op.opsm);
    glosize = [op.m ncols];

    spmd
        scheme  = sum(scheme);
        ypart   = ones(scheme,ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme;
        ycodist = codistributor1d(1,ygpart,glosize);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
    
else % Non-distributed
    y = ones(op.m,ncols);
end