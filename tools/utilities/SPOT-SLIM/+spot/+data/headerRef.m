function h = headerRef(header,index)
%HEADERREF  "Subsref" function for SeisDataContainer headers
%
%   h = headerRef(header,index) basically returns a header with a subsrefed
%   range (in the dimensional sense) of its header entities

% Index checking
index = index(1):index(end);
assert(isnumeric(index), 'index must be numeric!');
assert(isvector(index), 'index must be a vector!');
% singleton = false;

sizes   = header.size;
origin  = header.origin;
delta   = header.delta;
unit    = header.unit;
label   = header.label;

if length(origin) == 1 && all(index ~= 1) % indexing into singleton in vec case
    h        = header;
    h.dims   = 1;
    h.size   = 1;
    h.origin = 0;
    h.delta  = 1;
    h.unit   = {'u1'};
    h.label  = {'l1'};
else % not vec case
    % Check index range
    assert(index(1) >= 1 && index(end) <= length(sizes),...
        'index out of bounds');

    dims = length(sizes(index));

    % Fill in the new ranges
    h        = header;
    h.dims   = dims;
    h.size   = sizes(index);
    h.origin = origin(index);
    h.delta  = delta(index);
    h.unit   = unit(index);
    h.label  = label(index);
end % not vec case