function A = Src_Rcv_intep(mmz,mmx,mmy,z,x,y)
	% this function gives source and receiver interptation matrix "A"
	% apply A on source wavelet, it will put interptate off-grid source to the 
	% grid. A' will compute values off-grid with values on its surounding grid.
	
	

if nargin < 6, y = ones(size(x));end

np = length(z);

nz  = length(mmz);
nx  = length(mmx);
ny  = length(mmy);
mmz = mmz(:).';
mmx = mmx(:).';
mmy = mmy(:).';

n = ny * nz * nx;

if length(mmy) == 1
	z-x;% for check, z, x need to be same length
	A = spalloc(np,n,np*16);
	
	% dz = mmz(2) - mmz(1);
	% dx = mmx(2) - mmx(1);
	%
	% z_idx = floor((z - mmz(1))./dz);
	% z_temp = (z - mmz(z_idx))./dz;
	% ceoz   = lagrange_coe_special(z_temp);
	%
	%
	% x_idx  = floor((x - mmx(1))./dx);
	% x_temp = (x - mmx(x_idx))./dx;
	% ceox   = lagrange_coe_special(x_temp);
	%
	%
	% for mm = 1:np
	%
	%
	% end
	%

   

	for mm = 1:np
		
		dz = mmz(2) - mmz(1);
		dx = mmx(2) - mmx(1);
 	    z_id = floor((z(mm) - mmz(1))./dz)+1;
		x_id = floor((x(mm) - mmx(1))./dx)+1;
		
		switch x_id
			case 1
				x_id_grid = [x_id,x_id+1,x_id+2];
				ceos_x = lagrange_ceo(x(mm),[mmx(x_id)-dx,mmx(x_id_grid)]);
				ceos_x = ceos_x([2 3 4]);
			case nx-1
				x_id_grid = [x_id-1,x_id,x_id+1];
				ceos_x = lagrange_ceo(x(mm),[mmx(x_id_grid),mmx(end)+dx]);
				ceos_x = ceos_x([1 2 3]);
			case nx
				x_id_grid = [x_id-2,x_id-1,x_id];
				ceos_x = lagrange_ceo(x(mm),[mmx(x_id_grid),mmx(x_id)+dx]);
				ceos_x = ceos_x([1 2 3]);
			otherwise
				x_id_grid = [x_id-1,x_id,x_id+1,x_id+2];
				ceos_x = lagrange_ceo(x(mm),mmx(x_id_grid));
		end
		
		switch z_id
			case 1
				z_id_grid = [z_id,z_id+1,z_id+2];
				ceos_z = lagrange_ceo(z(mm),[mmz(z_id)-dz,mmz(z_id_grid)]);
				ceos_z = ceos_z([2 3 4]);
			case nz-1
				z_id_grid = [z_id-1,z_id,z_id+1];
				ceos_z = lagrange_ceo(z(mm),[mmz(z_id_grid),mmz(end)+dz]);
				ceos_z = ceos_z([1 2 3]);
			case nz
				z_id_grid = [z_id-2,z_id-1,z_id];
				ceos_z = lagrange_ceo(z(mm),[mmz(z_id_grid),mmz(z_id)+dz]);
				ceos_z = ceos_z([1 2 3]);
			otherwise
				z_id_grid = [z_id-1,z_id,z_id+1,z_id+2];
				ceos_z = lagrange_ceo(z(mm),mmz(z_id_grid));
		end
			
		ceos   = ceos_z*ceos_x';
		
		[XX,ZZ]   = meshgrid(x_id_grid,z_id_grid);
		ind		  = sub2ind([nz,nx,ny],ZZ,XX);

		A(mm,ind(:))  = ceos(:);
		

	end	
	
	
else
	z-x-y;% for check, z, x, y need to be same length
   A = spalloc(np,n,np*64);
   
   
	for mm = 1:np
	
		dz = mmz(2) - mmz(1);
		dx = mmx(2) - mmx(1);
		dy = mmy(2) - mmy(1);
	    z_id = floor((z(mm) - mmz(1))./dz)+1;
		x_id = floor((x(mm) - mmx(1))./dx)+1;
		y_id = floor((y(mm) - mmy(1))./dy)+1;
	
		
		switch x_id
			case 1
				x_id_grid = [x_id,x_id+1,x_id+2];
				ceos_x = lagrange_ceo(x(mm),[mmx(x_id)-dx,mmx(x_id_grid)]);
				ceos_x = ceos_x([2 3 4]);
			case nx-1
				x_id_grid = [x_id-1,x_id,x_id+1];
				ceos_x = lagrange_ceo(x(mm),[mmx(x_id_grid),mmx(end)+dx]);
				ceos_x = ceos_x([1 2 3]);
			case nx
				x_id_grid = [x_id-2,x_id-1,x_id];
				ceos_x = lagrange_ceo(x(mm),[mmx(x_id_grid),mmx(x_id)+dx]);
				ceos_x = ceos_x([1 2 3]);
			otherwise
				x_id_grid = [x_id-1,x_id,x_id+1,x_id+2];
				ceos_x = lagrange_ceo(x(mm),mmx(x_id_grid));
		end
	
		switch z_id
			case 1
				z_id_grid = [z_id,z_id+1,z_id+2];
				ceos_z = lagrange_ceo(z(mm),[mmz(z_id)-dz,mmz(z_id_grid)]);
				ceos_z = ceos_z([2 3 4]);
			case nz-1
				z_id_grid = [z_id-1,z_id,z_id+1];
				ceos_z = lagrange_ceo(z(mm),[mmz(z_id_grid),mmz(end)+dz]);
				ceos_z = ceos_z([1 2 3]);
			case nz
				z_id_grid = [z_id-2,z_id-1,z_id];
				ceos_z = lagrange_ceo(z(mm),[mmz(z_id_grid),mmz(z_id)+dz]);
				ceos_z = ceos_z([1 2 3]);
			otherwise
				z_id_grid = [z_id-1,z_id,z_id+1,z_id+2];
				ceos_z = lagrange_ceo(z(mm),mmz(z_id_grid));
		end
 
 
		switch y_id
			case 1
				y_id_grid = [y_id,y_id+1,y_id+2];
				ceos_y = lagrange_ceo(y(mm),[mmy(y_id)-dy,mmy(y_id_grid)]);
				ceos_y = ceos_y([2 3 4]);
			case ny-1
				y_id_grid = [y_id-1,y_id,y_id+1];
				ceos_y = lagrange_ceo(y(mm),[mmy(y_id_grid),mmy(end)+dy]);
				ceos_y = ceos_y([1 2 3]);
			case ny
				y_id_grid = [y_id-2,y_id-1,y_id];
				ceos_y = lagrange_ceo(y(mm),[mmy(y_id_grid),mmy(y_id)+dy]);
				ceos_y = ceos_y([1 2 3]);
			otherwise
				y_id_grid = [y_id-1,y_id,y_id+1,y_id+2];
				ceos_y = lagrange_ceo(y(mm),mmy(y_id_grid));
		end
		ceos   = kron(ceos_y',kron(ceos_x',ceos_z'));
     
		[XX,ZZ,YY]   = meshgrid(x_id_grid,z_id_grid,y_id_grid);
		ind		  = sub2ind([nz,nx,ny],ZZ,XX,YY);
	

		A(mm,ind(:)) = ceos(:);
   
	end
end