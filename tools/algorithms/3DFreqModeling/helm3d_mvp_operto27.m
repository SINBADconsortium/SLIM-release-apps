function y = helm3d_27pt_mvp(wn,h,n,npml,x,mode) %#codegen

MMM=1; NMM=2; PMM=3; MNM=4; NNM=5; PNM=6; MPM=7; NPM=8; PPM=9;
MMN=10;NMN=11;PMN=12;MNN=13;NNN=14;PNN=15;MPN=16;NPN=17;PPN=18;
MMP=19;NMP=20;PMP=21;MNP=22;NNP=23;PNP=24;MPP=25;NPP=26;PPP=27;

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
if mode==1
    for k=1:nz
        coef = complex(zeros(3,3,3),zeros(3,3,3));
        wn_window = complex(zeros(3,3,3),zeros(3,3,3));
        if k>1, if k < nz, zoffset = -1:1; else zoffset = -1:0; end, else zoffset = 0:1; end
        cxz = xz_coef*pmlz(k);
        cxzlo = xz_coef*pmlz_lo(k);
        cxzhi = xz_coef*pmlz_hi(k);
        cyz = yz_coef*pmlz(k);
        cyzlo = yz_coef*pmlz_lo(k);
        cyzhi = yz_coef*pmlz_hi(k);
        for j=1:ny
            if j > 1, if j < nz, yoffset = -1:1; else yoffset = -1:0; end, else yoffset = 0:1; end
            c1 = xy_coef*pmly(j) + cxz;
            c2_lo = cy*pmly_lo(j) + cyz;
            c2_hi = cy*pmly_hi(j) + cyz;
            c3_lo = yz_coef*pmly(j) + cz*pmlz_lo(k);
            c3_hi = yz_coef*pmly(j) + cz*pmlz_hi(k);
            for i=1:nx
                if i > 1, if i < nz, xoffset = -1:1; else xoffset = -1:0; end, else xoffset = 0:1; end
                % Load wn_window
                wn_window(2+xoffset,2+yoffset,2+zoffset) = wn(i+xoffset,j+yoffset,k+zoffset);
                
                % MNN
                coef(MNN) = cx*pmlx_lo(i) + c1 - wm2 * wn_window(1,2,2);
                % PNN
                coef(PNN) = cx*pmlx_hi(i) + c1 - wm2 * wn_window(3,2,2);
                % NMN
                coef(NMN) = xy_coef*pmlx(i) + c2_lo - wm2 * wn_window(2,1,2);
                % NPN
                coef(NPN) = xy_coef*pmlx(i) + c2_hi - wm2 * wn_window(2,3,2);
                % NNM
                coef(NNM) = xz_coef*pmlx(i) + c3_lo - wm2 * wn_window(2,2,1);
                % NNP
                coef(NNP) = xz_coef*pmlx(i) + c3_hi - wm2 * wn_window(2,2,3);
                % NMM
                coef(NMM) = 2*w3a*pmlx(i) - yz_coef*pmly_lo(j) - yz_coef*pmlz_lo(k) - wm3 * wn_window(2,1,1);
                % NMP
                coef(NMP) = 2*w3a*pmlx(i) - yz_coef*pmly_lo(j) - yz_coef*pmlz_hi(k) - wm3 * wn_window(2,1,3);
                % NPM
                coef(NPM) = 2*w3a*pmlx(i) - yz_coef*pmly_hi(j) - yz_coef*pmlz_lo(k) - wm3 * wn_window(2,3,1);
                % NPP
                coef(NPP) = 2*w3a*pmlx(i) - yz_coef*pmly_hi(j) - yz_coef*pmlz_hi(k) - wm3 * wn_window(2,3,3);
                % MNM
                coef(MNM) = -xz_coef*pmlx_lo(i) + 2*w3a*pmly(j) - xz_coef*pmlz_lo(k) - wm3 * wn_window(1,2,1);
                % MNP
                coef(MNP) = -xz_coef*pmlx_lo(i) + 2*w3a*pmly(j) - xz_coef*pmlz_hi(k) - wm3 * wn_window(1,2,3);
                % PNM
                coef(PNM) = -xz_coef*pmlx_hi(i) + 2*w3a*pmly(j) - xz_coef*pmlz_lo(k) - wm3 * wn_window(3,2,1);
                % PNP
                coef(PNP) = -xz_coef*pmlx_hi(i) + 2*w3a*pmly(j) - xz_coef*pmlz_hi(k) - wm3 * wn_window(3,2,3);
                % MMN
                coef(MMN) = -xy_coef*pmlx_lo(i) - xy_coef*pmly_lo(j) + 2*w3a*pmlz(k) - wm3 * wn_window(1,1,2);
                % MPN
                coef(MPN) = -xy_coef*pmlx_lo(i) - xy_coef*pmly_hi(j) + 2*w3a*pmlz(k) - wm3 * wn_window(1,3,2);
                % PMN
                coef(PMN) = -xy_coef*pmlx_hi(i) - xy_coef*pmly_lo(j) + 2*w3a*pmlz(k) - wm3 * wn_window(3,1,2);
                % PPN
                coef(PPN) = -xy_coef*pmlx_hi(i) - xy_coef*pmly_hi(j) + 2*w3a*pmlz(k) - wm3 * wn_window(3,3,2);
                % MMM
                coef(MMM) = -2*w3a*(pmlx_lo(i) + pmly_lo(j) + pmlz_lo(k)) - wm4 * wn_window(1,1,1);
                % MMP
                coef(MMP) = -2*w3a*(pmlx_lo(i) + pmly_lo(j) + pmlz_hi(k)) - wm4 * wn_window(1,1,3);
                % MPM
                coef(MPM) = -2*w3a*(pmlx_lo(i) + pmly_hi(j) + pmlz_lo(k)) - wm4 * wn_window(1,3,1);
                % MPP
                coef(MPP) = -2*w3a*(pmlx_lo(i) + pmly_hi(j) + pmlz_hi(k)) - wm4 * wn_window(1,3,3);
                % PMM
                coef(PMM) = -2*w3a*(pmlx_hi(i) + pmly_lo(j) + pmlz_lo(k)) - wm4 * wn_window(3,1,1);
                % PMP
                coef(PMP) = -2*w3a*(pmlx_hi(i) + pmly_lo(j) + pmlz_hi(k)) - wm4 * wn_window(3,1,3);
                % PPM
                coef(PPM) = -2*w3a*(pmlx_hi(i) + pmly_hi(j) + pmlz_lo(k)) - wm4 * wn_window(3,3,1);
                % PPP
                coef(PPP) = -2*w3a*(pmlx_hi(i) + pmly_hi(j) + pmlz_hi(k)) - wm4 * wn_window(3,3,3);
                % NNN
                coef(NNN) = -cx*pmlx(i) - cy*pmly(j) - cz*pmlz(k) + cNNN*wn_window(2,2,2);
                
                y(i,j,k) = sum(sum(sum( coef(2+xoffset,2+yoffset,2+zoffset) .* x(i+xoffset,j+yoffset,k+zoffset) )));
            end
        end
    end
    
end
y = vec(y);
end
