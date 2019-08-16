function update = Smoothupdate(up,mode)

mask = ones(size(up,1),1);
mask(1:16) = 0;
mask1 = smooth(smooth(mask));

if mode == 1
	% smooth mode
	update = repmat(mask1,1,size(up,2)).*up;
elseif mode == 2
	% hard mode
	update = repmat(mask,1,size(up,2)).*up;
end