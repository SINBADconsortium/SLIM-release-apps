function [outarr, mask] = hp1D_z (inarr, dz, dx, b_kz)
% high pass kz filtering
    
    [nz, nx] = size(inarr);
        
    MEarr = [inarr fliplr(inarr); flipud(inarr) rot90(inarr,2)];
    nz_ext = 2*nz;
    nx_ext = 2*nx;
    
    FTarr = fftshift(fft(MEarr));
    dk_x = 1/(dx*nx_ext);
    k_x = -dk_x*(ceil(nx_ext/2)-1):dk_x:dk_x*floor(nx_ext/2);
    dk_z = 1/(dz*nz_ext);
    k_z = -dk_z*(ceil(nz_ext/2)-1):dk_z:dk_z*floor(nz_ext/2);
    mask.k_x = k_x;
    mask.k_z = k_z;
    
    [x, z] = meshgrid(k_x, k_z);
    
    dist = abs(z./b_kz);
    
    mask.value = (cos((pi).*(dist - 1)./(1))+1)/2;
    
    mask.value(dist >= 2) = 0;
    mask.value(dist <= 1) = 1;
    mask.value = 1-mask.value;
    
    MEoutarr = ifft(ifftshift(mask.value .* FTarr));
    
    if isreal(inarr)
        outarr = real(double(MEoutarr(1:nz, 1:nx)));
        return
    end
    
    outarr = MEoutarr(1:nz,1:nx);

end