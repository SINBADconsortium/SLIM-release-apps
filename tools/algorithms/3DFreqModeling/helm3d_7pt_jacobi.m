function y = helm3d_7pt_jacobi(wn,h,n,npml,freq,alpha,b,x,mode) %#codegen

MMM=1; NMM=2; PMM=3; MNM=4; NNM=5; PNM=6; MPM=7; NPM=8; PPM=9;
MMN=10;NMN=11;PMN=12;MNN=13;NNN=14;PNN=15;MPN=16;NPN=17;PPN=18;
MMP=19;NMP=20;PMP=21;MNP=22;NNP=23;PNP=24;MPP=25;NPP=26;PPP=27;

nx = n(1); ny = n(2); nz = n(3);
hx2 = h(1)^2; hy2 = h(2)^2; hz2 = h(3)^2;

C = 10;
b = reshape(b,nx,ny,nz);
x = reshape(x,nx,ny,nz);
wn = reshape(wn,nx,ny,nz);

y = complex(zeros(nx,ny,nz),zeros(nx,ny,nz));

if mode==1
for k=1:nz
    coef = complex(zeros(3,3,3),zeros(3,3,3)); 
	if k>1, if k < nz, zoffset = -1:1; else zoffset = -1:0; end, else zoffset = 0:1; end
    pzM = pml_func1d(k-0.5,nz,npml(:,3),freq,C); pzN = pml_func1d(k,nz,npml(:,3),freq,C); pzP = pml_func1d(k+0.5,nz,npml(:,3),freq,C);
    pzlo = pzM*pzN/hz2; pzhi = pzN*pzP/hz2;
    coef(NNM) = -pzlo;
    coef(NNP) = -pzhi;
	for j=1:ny
		if j > 1, if j < nz, yoffset = -1:1; else yoffset = -1:0; end, else yoffset = 0:1; end		
        pyM = pml_func1d(j-0.5,ny,npml(:,2),freq,C); pyN = pml_func1d(j,ny,npml(:,2),freq,C); pyP = pml_func1d(j+0.5,ny,npml(:,2),freq,C);        
        pylo = pyM*pyN/hy2; pyhi = pyN*pyP/hy2;
        coef(NMN) = -pylo;
        coef(NPN) = -pyhi;
		for i=1:nx
            pxM = pml_func1d(i-0.5,nx,npml(:,1),freq,C); pxN = pml_func1d(i,nx,npml(:,1),freq,C); pxP = pml_func1d(i+0.5,nx,npml(:,1),freq,C);
            pxlo = pxM*pxN/hx2; pxhi = pxN*pxP/hx2;
            
			if i > 1, if i < nz, xoffset = -1:1; else xoffset = -1:0; end, else xoffset = 0:1; end
			
            coef(MNN) = -pxlo;
            coef(NNN) = pxlo + pxhi + pylo + pyhi + pzlo + pzhi - wn(i,j,k);
            coef(PNN) = -pxhi;

			Ax = sum(sum(sum( coef(2+xoffset,2+yoffset,2+zoffset) .* x(i+xoffset,j+yoffset,k+zoffset) )));
            r = b(i,j,k) - Ax;
            y(i,j,k) = x(i,j,k) + alpha*r/coef(NNN);
		end
	end
end

end
y = vec(y);
end

function func = pml_func1d(x,nx,nb,freq,C) %#codegen
    dist_from_int = (x <= nb(1)) .* (nb(1)-x) + (x > nx-nb(2)) .* (x-(nx-nb(2)+1));
    sigma = C * (dist_from_int/max(nb)).^2;
    func = (1-1i*sigma/real(freq)).^(-1);
end
