function v = fdct_v2c_mat(x,hdr)
%
%   Copyright 2008, Gilles Hennenfent

k = prod(hdr{1}{1});
v{1}{1} = reshape(x(1:k),hdr{1}{1});
for i=2:length(hdr)
   nw = length(hdr{i});

   for j = 1:nw
       ns = prod(hdr{i}{j});
       v{i}{j} = reshape(x(k+(1:ns)),hdr{i}{j});
       k = k + ns;
   end
end

