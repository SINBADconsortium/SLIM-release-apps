function y = rzeros(op,varargin)
%RZEROS Distributed zero vector in operator range
%
%   y = rzeros(A) generates a zero vector with the size of the operator
%   range, distributed accordingly to the operators needs so that A'*y is
%   a valid operation.
%
%   y = rzeros(A,NCOLS) generates NCOLS vectors

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
        ypart   = zeros(scheme,ncols);
        ygpart  = codistributed.zeros(1,numlabs);
        ygpart(labindex) = scheme;
        ycodist = codistributor1d(1,ygpart,glosize);
        y = codistributed.build(ypart,ycodist,'noCommunication');
    end
    
else % Non-distributed
    y = zeros(op.m,ncols);
end