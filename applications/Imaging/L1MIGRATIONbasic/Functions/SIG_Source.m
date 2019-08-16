function [b,nshot] = SIG_Source(xd,yd,nx,ny,mode,sx,sy,nx_border,ny_border);
% This function will generate sequential source function for YOGI and TIM's
% helmholz, if mode = 1; for YOGI. while if mode = 2 for TIM.
% 
% Input arguement:
%	xd, yd: mesh distant in lateral and depth dimension
%	nx, ny: grid number in lateral and depth dimension
%	sx, sy: shot position; if more than one shot
%			sx should be a vector
%	nx_border, ny_border: border grid numbers 
%	
if(nargin<10)
    nshot = 1;
end
if(nargin<9)
    ny_border = fix(ny*.1);
end
if(nargin<8)
    nx_border = fix(nx*.1);
end
if(nargin<7)
    sy = 5;
end
if(nargin<6)
    sx = round((nx+1)/2);
end

hx    = xd/(nx-1);
hy    = yd/(ny-1);
nshot = length(sx);

if mode == 1
    % b = sparse(zeros(nx*ny,1));   
    % ixys = (sy-1)*nx + sx;
    % b(ixys,1) = 1/(hx*hy);

	b = sparse(zeros(nx*ny,nshot));
	for m = 1:nshot
		ixys = (sy-1)*nx + sx(m);
		b(ixys,m) = 1/(hx*hy);
	end
elseif mode == 2
    % b = sparse(zeros((nx+2*nx_border)*(ny+2*ny_border),1));
    % ixys = (nx_border+sx-1)*(ny+2*ny_border)+ ny_border +sy;
    % b(ixys,1) = 1/(hx*hy);
	
	% b = sparse(zeros((nx+2*nx_border)*(ny+2*ny_border),nshot));
	b = spalloc((nx+2*nx_border)*(ny+2*ny_border),nshot,nshot);
	for m = 1:nshot
		ixys      = (nx_border+sx(m)-1)*(ny+2*ny_border)+ ny_border +sy;
		b(ixys,m) = 1/(hx*hy);
	end
end
