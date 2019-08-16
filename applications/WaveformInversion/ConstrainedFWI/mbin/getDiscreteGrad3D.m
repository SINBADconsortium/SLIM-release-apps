function [Dx_3D,Dy_3D,Dz_3D] = getDiscreteGrad3D(n,d)
%depth dimension (z) is the 3rd dimension
n1=n(1);n2=n(2);n3=n(3);
h1=d(1);h2=d(2);h3=d(3);


%define difference matrix D acting on vectorized model using kron
I1 = speye(n1); %x
I2 = speye(n2); %y
I3 = speye(n3); %z
Dz = spdiags(ones(n3,1)*[1 -1],-1:0,n3,n3-1)/h3;
Dx = spdiags(ones(n1,1)*[1 -1],-1:0,n1,n1-1)/h1;
Dy = spdiags(ones(n2,1)*[1 -1],-1:0,n2,n2-1)/h2;

%D = [kron(kron(I3,I2),Dz) ;
%     kron(kron(I3,Dz),I1) ;
%     kron(Dz,kron(I2,I1))] ;
 
Dx_3D=kron(kron(I3,I2),Dx)';
Dy_3D= kron(kron(I3,Dy),I1)';
Dz_3D= kron(Dz,kron(I2,I1))';
     
end