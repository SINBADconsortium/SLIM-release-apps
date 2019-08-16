function y = Helm3dmvp_v(wn,h,n,npml,x,mode) %#codegen

nx = n(1); ny = n(2); nz = n(3);
hx2 = h(1)^2; hy2 = h(2)^2; hz2 = h(3)^2;
hyz  = hy2 + hz2; hxy  = hy2 + hx2; hxz  = hx2 + hz2; hxyz = hx2 + hy2 + hz2;

x = reshape(x,nx,ny,nz);
wn = reshape(wn,nx,ny,nz);

% pml handling
[pmlx, pmly, pmlz] = create_pml(nx, ny, nz, npml);
pmlx_lo = pmlx(:,1); pmlx_hi = pmlx(:,2); pmlx = pmlx_lo + pmlx_hi;
pmly_lo = pmly(:,1); pmly_hi = pmly(:,2); pmly = pmly_lo + pmly_hi;
pmlz_lo = pmlz(:,1); pmlz_hi = pmlz(:,2); pmlz = pmlz_lo + pmlz_hi;

% weights
w1  = 1.8395262e-5;
w2  = (0.890077)/3;
w3  = (0.1099046)/4;
wm1 = 0.49649658;
wm2 = (0.4510125)/6;
wm3 = (0.052487)/12;
wm4 = (0.45523e-5)/8;
w3a = ((w3*3)/(4*hxyz));

% static coefficients
cx = - (w1/hx2 + w2/hx2 + w2/hxz + w2/hxy + 8*w3a);
cy = - (w1/hy2 + w2/hy2 + w2/hyz + w2/hxy + 8*w3a);
cz = - (w1/hz2 + w2/hz2 + w2/hxz + w2/hyz + 8*w3a);
cNNN = - (w1 + 3*w2 + (16*w3a*hxyz)/3 + wm1-1);
y = complex(zeros(nx,ny,nz),zeros(nx,ny,nz));

xz_coef = w2/(2*hxz);
xy_coef = w2/(2*hxy);
yz_coef = w2/(2*hyz);
xN = 1:nx; xP = 2:nx; xM = 1:nx-1;
yN = 1:ny; yP = 2:ny; yM = 1:ny-1;
zN = 1:nz; zP = 2:nz; zM = 1:nz-1;
if mode==1
	% MNN
	coef = xyz_sum(cx*pmlx_lo , xy_coef*pmly , xz_coef*pmlz);
	coef(xP,yN,zN) = coef(xP,yN,zN) - wm2*wn(xM,yN,zN);
	y(xP,yN,zN) = y(xP,yN,zN) + coef(xP,yN,zN) .* x(xM,yN,zN);

	% PNN
	coef = xyz_sum(cx*pmlx_hi , xy_coef*pmly , xz_coef*pmlz);
	coef(xM,yN,zN) = coef(xM,yN,zN) - wm2*wn(xP,yN,zN);
	y(xM,yN,zN) = y(xM,yN,zN) + coef(xM,yN,zN) .* x(xP,yN,zN);

	% NMN
	coef = xyz_sum(xy_coef*pmlx , cy*pmly_lo , yz_coef*pmlz);
	coef(xN,yP,zN) = coef(xN,yP,zN) - wm2*wn(xN,yM,zN);
	y(xN,yP,zN) = y(xN,yP,zN) + coef(xN,yP,zN) .* x(xN,yM,zN);

	% NPN
	coef = xyz_sum(xy_coef*pmlx , cy*pmly_hi , yz_coef*pmlz);
	coef(xN,yM,zN) = coef(xN,yM,zN) - wm2*wn(xN,yP,zN);
	y(xN,yM,zN) = y(xN,yM,zN) + coef(xN,yM,zN) .* x(xN,yP,zN);

	% NNM
	coef = xyz_sum(xz_coef*pmlx , yz_coef*pmly , cz*pmlz_lo);
	coef(xN,yN,zP) = coef(xN,yN,zP) - wm2*wn(xN,yN,zM);
	y(xN,yN,zP) = y(xN,yN,zP) + coef(xN,yN,zP) .* x(xN,yN,zM);

	% NNP
	coef = xyz_sum(xz_coef*pmlx , yz_coef*pmly , cz*pmlz_hi);
	coef(xN,yN,zM) = coef(xN,yN,zM) - wm2*wn(xN,yN,zP);
	y(xN,yN,zM) = y(xN,yN,zM) + coef(xN,yN,zM) .* x(xN,yN,zP);

	% NMM
	coef = xyz_sum(2*w3a*pmlx , -yz_coef*pmly_lo , -yz_coef*pmlz_lo);
	coef(xN,yP,zP) = coef(xN,yP,zP) - wm3*wn(xN,yM,zM);
	y(xN,yP,zP) = y(xN,yP,zP) + coef(xN,yP,zP) .* x(xN,yM,zM);

	% NMP
	coef = xyz_sum(2*w3a*pmlx , -yz_coef*pmly_lo , -yz_coef*pmlz_hi);
	coef(xN,yP,zM) = coef(xN,yP,zM) - wm3*wn(xN,yM,zP);
	y(xN,yP,zM) = y(xN,yP,zM) + coef(xN,yP,zM) .* x(xN,yM,zP);

	% NPM
	coef = xyz_sum(2*w3a*pmlx , -yz_coef*pmly_hi , -yz_coef*pmlz_lo);
	coef(xN,yM,zP) = coef(xN,yM,zP) - wm3*wn(xN,yP,zM);
	y(xN,yM,zP) = y(xN,yM,zP) + coef(xN,yM,zP) .* x(xN,yP,zM);

	% NPP
	coef = xyz_sum(2*w3a*pmlx , -yz_coef*pmly_hi , -yz_coef*pmlz_hi);
	coef(xN,yM,zM) = coef(xN,yM,zM) - wm3*wn(xN,yP,zP);
	y(xN,yM,zM) = y(xN,yM,zM) + coef(xN,yM,zM) .* x(xN,yP,zP);

	% MNM
	coef = xyz_sum(-xz_coef*pmlx_lo , 2*w3a*pmly , -xz_coef*pmlz_lo);
	coef(xP,yN,zP) = coef(xP,yN,zP) - wm3*wn(xM,yN,zM);
	y(xP,yN,zP) = y(xP,yN,zP) + coef(xP,yN,zP) .* x(xM,yN,zM);

	% MNP
	coef = xyz_sum(-xz_coef*pmlx_lo , 2*w3a*pmly , -xz_coef*pmlz_hi);
	coef(xP,yN,zM) = coef(xP,yN,zM) - wm3*wn(xM,yN,zP);
	y(xP,yN,zM) = y(xP,yN,zM) + coef(xP,yN,zM) .* x(xM,yN,zP);

	% PNM
	coef = xyz_sum(-xz_coef*pmlx_hi , 2*w3a*pmly , -xz_coef*pmlz_lo);
	coef(xM,yN,zP) = coef(xM,yN,zP) - wm3*wn(xP,yN,zM);
	y(xM,yN,zP) = y(xM,yN,zP) + coef(xM,yN,zP) .* x(xP,yN,zM);

	% PNP
	coef = xyz_sum(-xz_coef*pmlx_hi , 2*w3a*pmly , -xz_coef*pmlz_hi);
	coef(xM,yN,zM) = coef(xM,yN,zM) - wm3*wn(xP,yN,zP);
	y(xM,yN,zM) = y(xM,yN,zM) + coef(xM,yN,zM) .* x(xP,yN,zP);

	% MMN
	coef = xyz_sum(-xy_coef*pmlx_lo , -xy_coef*pmly_lo , 2*w3a*pmlz);
	coef(xP,yP,zN) = coef(xP,yP,zN) - wm3*wn(xM,yM,zN);
	y(xP,yP,zN) = y(xP,yP,zN) + coef(xP,yP,zN) .* x(xM,yM,zN);

	% MPN
	coef = xyz_sum(-xy_coef*pmlx_lo , -xy_coef*pmly_hi , 2*w3a*pmlz);
	coef(xP,yM,zN) = coef(xP,yM,zN) - wm3*wn(xM,yP,zN);
	y(xP,yM,zN) = y(xP,yM,zN) + coef(xP,yM,zN) .* x(xM,yP,zN);

	% PMN
	coef = xyz_sum(-xy_coef*pmlx_hi , -xy_coef*pmly_lo , 2*w3a*pmlz);
	coef(xM,yP,zN) = coef(xM,yP,zN) - wm3*wn(xP,yM,zN);
	y(xM,yP,zN) = y(xM,yP,zN) + coef(xM,yP,zN) .* x(xP,yM,zN);

	% PPN
	coef = xyz_sum(-xy_coef*pmlx_hi , -xy_coef*pmly_hi , 2*w3a*pmlz);
	coef(xM,yM,zN) = coef(xM,yM,zN) - wm3*wn(xP,yP,zN);
	y(xM,yM,zN) = y(xM,yM,zN) + coef(xM,yM,zN) .* x(xP,yP,zN);

	% MMM
	coef = xyz_sum(-2*w3a*pmlx_lo , -2*w3a*pmly_lo , -2*w3a*pmlz_lo);
	coef(xP,yP,zP) = coef(xP,yP,zP) - wm4*wn(xM,yM,zM);
	y(xP,yP,zP) = y(xP,yP,zP) + coef(xP,yP,zP) .* x(xM,yM,zM);

	% MMP
	coef = xyz_sum(-2*w3a*pmlx_lo , -2*w3a*pmly_lo , -2*w3a*pmlz_hi);
	coef(xP,yP,zM) = coef(xP,yP,zM) - wm4*wn(xM,yM,zP);
	y(xP,yP,zM) = y(xP,yP,zM) + coef(xP,yP,zM) .* x(xM,yM,zP);

	% MPM
	coef = xyz_sum(-2*w3a*pmlx_lo , -2*w3a*pmly_hi , -2*w3a*pmlz_lo);
	coef(xP,yM,zP) = coef(xP,yM,zP) - wm4*wn(xM,yP,zM);
	y(xP,yM,zP) = y(xP,yM,zP) + coef(xP,yM,zP) .* x(xM,yP,zM);

	% MPP
	coef = xyz_sum(-2*w3a*pmlx_lo , -2*w3a*pmly_hi , -2*w3a*pmlz_hi);
	coef(xP,yM,zM) = coef(xP,yM,zM) - wm4*wn(xM,yP,zP);
	y(xP,yM,zM) = y(xP,yM,zM) + coef(xP,yM,zM) .* x(xM,yP,zP);

	% PMM
	coef = xyz_sum(-2*w3a*pmlx_hi , -2*w3a*pmly_lo , -2*w3a*pmlz_lo);
	coef(xM,yP,zP) = coef(xM,yP,zP) - wm4*wn(xP,yM,zM);
	y(xM,yP,zP) = y(xM,yP,zP) + coef(xM,yP,zP) .* x(xP,yM,zM);

	% PMP
	coef = xyz_sum(-2*w3a*pmlx_hi , -2*w3a*pmly_lo , -2*w3a*pmlz_hi);
	coef(xM,yP,zM) = coef(xM,yP,zM) - wm4*wn(xP,yM,zP);
	y(xM,yP,zM) = y(xM,yP,zM) + coef(xM,yP,zM) .* x(xP,yM,zP);

	% PPM
	coef = xyz_sum(-2*w3a*pmlx_hi , -2*w3a*pmly_hi , -2*w3a*pmlz_lo);
	coef(xM,yM,zP) = coef(xM,yM,zP) - wm4*wn(xP,yP,zM);
	y(xM,yM,zP) = y(xM,yM,zP) + coef(xM,yM,zP) .* x(xP,yP,zM);

	% PPP
	coef = xyz_sum(-2*w3a*pmlx_hi , -2*w3a*pmly_hi , -2*w3a*pmlz_hi);
	coef(xM,yM,zM) = coef(xM,yM,zM) - wm4*wn(xP,yP,zP);
	y(xM,yM,zM) = y(xM,yM,zM) + coef(xM,yM,zM) .* x(xP,yP,zP);

	% NNN
	coef = xyz_sum(-cx*pmlx , -cy*pmly , -cz*pmlz);
	coef = coef + cNNN*wn(xN,yN,zN);
	y(xN,yN,zN) = y(xN,yN,zN) + coef(xN,yN,zN) .* x(xN,yN,zN);
end