function M = src_mask( model, r )
%SRCREC_MASK Summary of this function goes here
%   Detailed explanation goes here

M = ones(model.n);
nsx = length(model.xsrc); nsy = length(model.ysrc); nsz = length(model.zsrc);
if nsx == 1, fixed_src_dim = 1; elseif nsy==1, fixed_src_dim = 2; elseif nsz==1, fixed_src_dim = 3; end
nsrc = nsx*nsy*nsz;
[x,y,z] = odn2grid(model.o,model.d,model.n);
[X,Y,Z] = ndgrid(x,y,z);
for i=1:nsrc   
   switch fixed_src_dim
       case 1
           srcx = 1;
           [srcy,srcz] = ind2sub([nsy,nsz],i);
       case 2
           srcy = 1;
           [srcx,srcz] = ind2sub([nsx,nsz],i);
       case 3
           srcz = 1;           
           [srcx,srcy] = ind2sub([nsx,nsy],i);
   end
   M((X-model.xsrc(srcx)).^2 + (Y-model.ysrc(srcy)).^2 + (Z-model.zsrc(srcz)).^2 <= r^2 ) = 0;
end

end

