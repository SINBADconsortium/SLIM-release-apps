function h = headerMod(op,header,mode)
%HEADERMOD  Header Modifying Function for SeisDataContainer compatibility
%
%   h = headerMod(op,xmeta,header,mode) returns a modified header 
%   corresponding to the actions the spot operator would have on the 
%   metadata of the SeisDataContainer.
%
%   op = Spot operator
%   xmeta = explicit metadata stored on the datacon
%   header = header of the

% By default this only changes the implicit header size

% Extract explicit size indices
exsize = header.exsize;

if mode == 1
    h = header; % Copy header
    % Replace old first (collapsed) dimensional sizes with operator sizes.
    h.size(exsize(1,1):exsize(2,1)) = [];
    h.size = [op.ms{:} h.size];
else
    h = header;
    h.size(exsize(1,1):exsize(2,1)) = [];
    h.size = [op.ns{:} h.size];
end

exsize_out = 1:length(h.size);
exsize_out = [exsize_out;exsize_out];
h.exsize   = exsize_out;