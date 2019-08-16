function Data = DistPermute(Data, perm, dist, varargin)
% DISTPERMUTE(Data, perm, dist) - to permute distributed nd-arrays. Data is
% a distributed array, perm is the desired arrangement of dimensions, and
% the optional dist specifies the dimension from the original array for the
% output to be distributed along. If dist isn't provided, the output is
% distributed along its last dimension.

numdims = ndims(Data);
if ~(nargin == 5 && ischar( varargin{1} ) && ...
        strcmp( varargin{1}, 'internalcall' ))
    %no need to check inputs for internal calls
    narginchk(2,3)
    if nargin == 2,     dist = perm(numdims);  end
    if numlabs ~= 1
        error('This function is not meant to be called from an spmd block')
    end
    if ~isa(Data, 'distributed')
        error('Input must be distributed')
    end
    if ~isnumeric(dist) || length(dist) ~= 1 || dist > numdims
        error('Distribution dimension must be an individual dimension of Data')
    end
    temp = perms(1:numdims);
    for i = 1:length(temp)
        if perm == temp(i,:), break;  end
        if i == length(temp),   error('Invalid permutation vector'),   end
    end
    
    %get dimensions of the data and rearrange to the dimensions after the permute
    temp = size(Data);
    dims = zeros(1,numdims);
    for i = 1:numdims
        dims(i) = temp(perm(i));
    end
else
    dims = varargin{2};
end

spmd
    codistr = getCodistributor(Data);
    dimdist = codistr.Dimension;
    
    %incase the distributed dimension is involved in the permute
    dimdist = find(perm == dimdist);
    dist = find(perm == dist);
    
    %do the permute locally on each lab and then build the codistributed
    %array again
    Data = permute(getLocalPart(Data), perm);
    codistr = codistributor1d( dimdist, ...
        codistributor1d.unsetPartition, dims);
    Data = codistributed.build(Data, codistr);
    
    %if the dimension of distribution changed
    if dimdist ~= dist
        codistr = codistributor1d(dist);
        Data = redistribute( Data, codistr);
    end
end
end
