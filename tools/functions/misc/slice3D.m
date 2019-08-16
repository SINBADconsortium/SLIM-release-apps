function slice3D( D, grid_struct,ix, iy, iz )
%SLICE3D Display slices of a 3D volume.
%  
%  Curt Da Silva, 2015
%
%  Usage:
%    slice3D( D, grid_struct, ix, iy, iz );
% 
%  Input:
%    D           - 3D volume conforming to grid_struct parameters
%    grid_struct - struct with .{o,d,n} fields that determine
%                  {x,y,z} coordinates of the grid
%    ix, iy, iz  - coordinates corresponding to the planes to display
%    
%  Example:
%    slice3D( D, grid_struct, 100, [],[] ); - plots a x=100 slice
%    slice3D( D, grid_struct, [100 500], [], []); plots a x=100 and x=500 slice
%    slice3D( D, grid_struct, 100, [], 900 ); - plots a x=100 slice and a z=900 slice
   


if isfield(grid_struct,'nt')
    o = grid_struct.ot; d = grid_struct.dt; n = grid_struct.nt;
elseif isfield(grid_struct,'n')
    o = grid_struct.o; d = grid_struct.d; n = grid_struct.n;
else
    error('Invalid grid struct');
end

D = reshape(D,n);
[x,y,z] = odn2grid(o,d,n);
z = flip(z);
D = flip(D,3);
[X,Y,Z] = meshgrid(x,y,z);
D = permute(D,[2 1 3]);

figure;
slice(X,Y,Z,D,ix,iy,iz);
set(gca,'ZDir','reverse');

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
% removes grid lines
set(findobj(gca,'Type','Surface'),'EdgeColor','none')
set(gcf,'Color','w');
end

