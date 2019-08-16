function y = helm3d_mvp_chen27(v,freq,h,n,npml,x,mode) %#codegen

MMM=1; NMM=2; PMM=3; MNM=4; NNM=5; PNM=6; MPM=7; NPM=8; PPM=9;
MMN=10;NMN=11;PMN=12;MNN=13;NNN=14;PNN=15;MPN=16;NPN=17;PPN=18;
MMP=19;NMP=20;PMP=21;MNP=22;NNP=23;PNP=24;MPP=25;NPP=26;PPP=27;

nx = n(1); ny = n(2); nz = n(3);
dx = h(1)^2; dy = h(2)^2; dz = h(3)^2;
cx = 1/dx; cy = 1/dy; cz = 1/dz;
omega2 = (2*pi*freq)^2;

[~,I] = max(h);
switch I
  case 1
    r1 = h(1)/h(2);
    r2 = h(1)/h(3);
  case 2
    r1 = h(2)/h(1);
    r2 = h(2)/h(3);
  case 3
    r1 = h(3)/h(1);
    r2 = h(3)/h(2);
  otherwise
    error('Error');
end


if r1==1
    if r2==1
        alpha1 = 0.097426; alpha2 = 0.001449; beta1 = 0.042419; beta2 = 0.028953; gamma1 = 0.100520; gamma2 = 0.000000; c = 0.474309; d = 0.084057; e = 0.001779;    
    elseif r2==2
        alpha1 = 0.027756; alpha2 = 0.035823; beta1 = 0.099402; beta2 = 0.000000; gamma1 = 0.075170; gamma2 = 0.010711; c = 0.468043; d = 0.086909; e = 0.000875; 
    elseif r2==3
        alpha1 = 0.091501; alpha2 = 0.005050; beta1 = 0.101582; beta2 = 0.000000; gamma1 = 0.062715; gamma2 = 0.016398; c = 0.454915; d = 0.090845; e = 0.000000;
    else
        error('Unsupported aspect ratio');
    end
else
    error('Error');
end
%alpha1 = 0; alpha2 = 0; beta1 = 0; beta2 = 0; gamma1 = 0; gamma2 = 0; c = 1; d = 0; e = 0;
f = (1-c-6*d-12*e)/8;

alpha1 = alpha1/dx;
alpha2 = alpha2/dx;
beta1 = beta1/dy;
beta2 = beta2/dy;
gamma1 = gamma1/dz;
gamma2 = gamma2/dz;

v = reshape(v,nx,ny,nz);
x = reshape(x,nx,ny,nz);

pmlx = pml_func1d(1:nx,nx,npml(:,1),freq);
pmlxlo = pml_func1d((1:nx)-0.5,nx,npml(:,1),freq);
pmlxhi = pml_func1d((1:nx)+0.5,nx,npml(:,1),freq);

pmly = pml_func1d(1:ny,ny,npml(:,2),freq);
pmlylo = pml_func1d((1:ny)-0.5,ny,npml(:,2),freq);
pmlyhi = pml_func1d((1:ny)+0.5,ny,npml(:,2),freq);

pmlz = pml_func1d(1:nz,nz,npml(:,1),freq);
pmlzlo = pml_func1d((1:nz)-0.5,nz,npml(:,1),freq);
pmlzhi = pml_func1d((1:nz)+0.5,nz,npml(:,1),freq);

y = complex(zeros(nx,ny,nz),zeros(nx,ny,nz));

if mode==1
    parfor k=1:nz
        coef = complex(zeros(3,3,3),zeros(3,3,3));
        if k>1, if k < nz, zoffset = -1:1; else zoffset = -1:0; end, else zoffset = 0:1; end
        pzlo = pmlz(k)*pmlzlo(k); pzhi = pmlz(k)*pmlzhi(k);    
        g1_lo = gamma1*pzlo; g2_lo = gamma2*pzlo;
        g1_hi = gamma1*pzhi; g2_hi = gamma2*pzhi;
        for j=1:ny
            if j > 1, if j < nz, yoffset = -1:1; else yoffset = -1:0; end, else yoffset = 0:1; end
            pylo = pmly(j)*pmlylo(j); pyhi = pmly(j)*pmlyhi(j);    
            b1_lo = beta1*pylo; b2_lo = beta2*pylo;
            b1_hi = beta1*pyhi; b2_hi = beta2*pyhi;            
            for i=1:nx
                if i > 1, if i < nz, xoffset = -1:1; else xoffset = -1:0; end, else xoffset = 0:1; end
                pxlo = pmlx(i)*pmlxlo(i); pxhi = pmlx(i)*pmlxhi(i);    
                a1_lo = alpha1*pxlo; a2_lo = alpha2*pxlo;
                a1_hi = alpha1*pxhi; a2_hi = alpha1*pxhi;
                
                % Load wn_window
                wn = omega2/v(i,j,k)^2;
                
                fwn = f*wn; ewn = e*wn;
                
                coef(MMM) = a2_lo + b2_lo + fwn + g2_lo;
                coef(MMN) = a1_lo + b1_lo + ewn - g2_hi - g2_lo;
                coef(MMP) = a2_lo + b2_lo + fwn + g2_hi;
                coef(MNM) = a1_lo - b2_hi - b2_lo + ewn + g1_lo;
                coef(MNN) = -4*a1_lo - 4*a2_lo - b1_hi - b1_lo + d*wn - g1_hi - g1_lo + pxlo*cx;
                coef(MNP) = a1_lo - b2_hi - b2_lo + ewn + g1_hi;
                coef(MPM) = a2_lo + b2_hi + fwn + g2_lo;
                coef(MPN) = a1_lo + b1_hi + ewn - g2_hi - g2_lo;
                coef(MPP) = a2_lo + b2_hi + fwn + g2_hi;
                coef(NMM) = -a2_hi - a2_lo + b1_lo + ewn + g1_lo;
                coef(NMN) = -a1_hi - a1_lo - 4*b1_lo - 4*b2_lo + d*wn - g1_hi - g1_lo + pylo*cy;
                coef(NMP) = -a2_hi - a2_lo + b1_lo + ewn + g1_hi;
                coef(NNM) = -a1_hi - a1_lo - b1_hi - b1_lo + d*wn - 4*g1_lo - 4*g2_lo + pzlo*cz;
                coef(NNN) = 4*a1_hi + 4*a1_lo + 4*a2_hi + 4*a2_lo + 4*b1_hi + 4*b1_lo + 4*b2_hi + 4*b2_lo + c*wn + 4*g1_hi + 4*g1_lo + 4*g2_hi + 4*g2_lo - pzhi*cz - pzlo*cz - pyhi*cy - pylo*cy - pxhi*cx - pxlo*cx;
                coef(NNP) = -a1_hi - a1_lo - b1_hi - b1_lo + d*wn - 4*g1_hi - 4*g2_hi + pzhi*cz;
                coef(NPM) = -a2_hi - a2_lo + b1_hi + ewn + g1_lo;
                coef(NPN) = -a1_hi - a1_lo - 4*b1_hi - 4*b2_hi + d*wn - g1_hi - g1_lo + pyhi*cy;
                coef(NPP) = -a2_hi - a2_lo + b1_hi + ewn + g1_hi;
                coef(PMM) = a2_hi + b2_lo + fwn + g2_lo;
                coef(PMN) = a1_hi + b1_lo + ewn - g2_hi - g2_lo;
                coef(PMP) = a2_hi + b2_lo + fwn + g2_hi;
                coef(PNM) = a1_hi - b2_hi - b2_lo + ewn + g1_lo;
                coef(PNN) = -4*a1_hi - 4*a2_hi - b1_hi - b1_lo + d*wn - g1_hi - g1_lo + pxhi*cx;
                coef(PNP) = a1_hi - b2_hi - b2_lo + ewn + g1_hi;
                coef(PPM) = a2_hi + b2_hi + fwn + g2_lo;
                coef(PPN) = a1_hi + b1_hi + ewn - g2_hi - g2_lo;
                coef(PPP) = a2_hi + b2_hi + fwn + g2_hi;                                
                
                y(i,j,k) = sum(sum(sum( (coef(2+xoffset,2+yoffset,2+zoffset)) .* x(i+xoffset,j+yoffset,k+zoffset) )));
            end
        end
    end
    
end
y = vec(y);
end


function func = pml_func1d(x,nx,nb,freq) %#codegen
    C = 10;
    dist_from_int = (x <= nb(1)) .* (nb(1)-x) + (x > nx-nb(2)) .* (x-(nx-nb(2)+1));
    sigma = C * (dist_from_int/max(nb)).^2;
    func = (1-1i*sigma/freq).^(-1);
end
