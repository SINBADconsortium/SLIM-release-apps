function Mask = make_mask(Nt,Nr,top,t_offset,slope)
% function that creates a triangular mask, written for opTriMute

% make the logical mask matrix
Mask = repmat(logical(0),Nt,Nr);

% top
Mask(1:top,:) = 1;

% 
for x = 1:Nr
    t = ceil(t_offset + slope * (x-1)) - 1;
    Mask(([1:Nt] < t),x) = 1;
end

Mask2 = fliplr(Mask);

Mask = [Mask2 Mask(:,[2:end])];

end