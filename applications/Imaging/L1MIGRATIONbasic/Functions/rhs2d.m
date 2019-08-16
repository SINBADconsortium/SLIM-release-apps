function [b,ixyshot] = rhs2d(xd,yd,nx,ny,nshot)

b = sparse(zeros(nx*ny,nshot));

hx = xd/(nx+1);
hy = yd/(ny+1);

if nshot == 1
  ixs = round((nx+1)/2);
  iys = 3;
  ixys = (iys-1)*nx + ixs;
  b(ixys,nshot) = 1/(hx*hy);
  ixyshot(1) = ixs;
end

%if nshot == 2
%  ixs1 = 1+3;
%  iys1 = 3;
%  ixys1 = (iys1-1)*nx + ixs1;
%  b(ixys1,1) = 1/(hx*hy);
%  ixyshot(1) = ixs1;
%  
%  ixs2 = nx-3;
%  iys2 = 3;
%  ixys2 = (iys2-1)*nx + ixs2;
%  b(ixys2,2) = 1/(hx*hy);
%  ixyshot(2) = ixs2;
%end

if nshot == 4
  for i = 1:nshot
     ixs = (i-1)*74+74;
     iys = 3;
     ixys = (iys-1)*nx+ixs;
     b(ixys,i) = 1/(hx*hy);
     ixyshot(i) = ixys;
  end
end

if nshot == 37
  for i = 1:nshot
    ixs = (i-1)*10+5;
    iys = 3;
    ixys = (iys-1)*nx+ixs;
    b(ixys,1) = 1;
    ixyshot(i) = ixys;
  end
end
    
if nshot == round((nx+1)/2)
  for i = 1:nshot
    ixs = (i-1)*2+1;
    iys = 3;
    ixys = (iys-1)*nx + ixs;
    b(ixys,i) = 1/(hx*hy);
    ixyshot(i) = ixs;
  end
end

if nshot == nx
  for ixs = 1:nx
    iys = 3;
    ixys = (iys-1)*nx + ixs;
    b(ixys,ixs) = 1/(hx*hy);
    ixyshot(ixs) = ixs;
  end
end


whos b

