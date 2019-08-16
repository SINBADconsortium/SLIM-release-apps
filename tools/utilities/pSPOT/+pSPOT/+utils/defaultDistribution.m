function part = defaultDistribution(dsize)
%DEFAULTDISTRIBUTION Generates the default partition size for your
%                   distribtued array
%
%   PARTITION = defaultDistribution(DSIZE) generates a vector containing
%   the default partition based on the size of the desired distribution 
%   dimension.
%
%   Note that this function must be called from without a spmd block.
%   
%   See codistributor1d.defaultPartition for more details

nlabs = parpool_size();
r     = rem(dsize,nlabs);
c     = ceil(dsize/nlabs);
f     = floor(dsize/nlabs);
part(1:r)       = c;
part(r+1:nlabs) = f;
