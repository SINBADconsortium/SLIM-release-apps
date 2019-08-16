function [H,dH,ddH,A] = Helm2D_opt(k,d,n,nb,unit,freq,f0,diag_offset)
% HELM2D_OPT - Optimized 9pt scheme for 2D Helmholtz, from 
%    'AN OPTIMAL 9-POINT FINITE DIFFERENCE SCHEME FOR THE HELMHOLTZ
%    EQUATION WITH PML', Z. CHEN, et al. 2013
%
%
% Curt Da Silva, 2016 
%
% Usage:
%   [H,dH,ddH,A] = Helm2D_opt(k,d,n,nb,unit,freq,f0,diag_offset);
%
% Input:
%   k     - gridded wavespeed or slowness^2 parameter with pml
%   d     - grid spacing
%   n     - model size (with pml)
%   nb    - number of pml points in each dimension (2 x 3 matrix)
%   freq  - Helmholtz frequency
%   f0    - peak frequency of source wavelet (used in pml)
%
% Output:
%   H     - sparse Helmholtz matrix
%  dH     - derivative of Helmholtz wrto model parameter
%  ddH    - second derivative of Helmholtz wrto model parameter
%   A     - 

    if exist('diag_offset','var')==0, diag_offset = 0; end
    assert(size(nb,1)==2 && size(nb,2)==2,'nb must be a 2 x 2 matrix');
% Convert model parameter to rad^2*(s^2/m^2)    
    [fm,dfm,ddfm] = param2wavenum(k,freq,unit);
    fm = fm + diag_offset;
    nz = n(1); nx = n(2); N = nz*nx;
    hz = d(1); hx = d(2);
    a0 = 1.79;
    
    % Choose optimal averaging constants for laplacian, averaging matrices
    b = 0.7926;
    d = 0.3768; e = -0.0064; 
    
    consts = [0.6803,0.4444,0.0008;
              0.7427,0.4088,-0.0036;
              0.7840,0.3832,-0.0060;
              0.8020,0.3712,-0.0072;
              0.8133,0.3637,-0.0075;
              0.8219,0.3578,-0.0078;
              0.8271,0.3540,-0.0080];
    if freq >= 2.5 && freq < 3, i = 1; 
    elseif freq >= 3 && freq < 4, i = 2;
    elseif freq >= 4 && freq < 5, i = 3;
    elseif freq >= 5 && freq < 6, i = 4;
    elseif freq >= 6 && freq < 8, i = 5;
    elseif freq >= 8 && freq < 10, i = 6;
    elseif freq > 10, i = 7; 
    else i = [];
    end
    
    if ~isempty(i), b = consts(i,1); d = consts(i,2); e = consts(i,3); end
    %b = 1; d = 0; e = 0;
    c = 1-d-e;    
    
    sdiag = @(x) spdiags(vec(x),0,N,N);
    
    ez = pml_func1d(nz,nb(:,1),a0,f0,freq);
    ex = pml_func1d(nx,nb(:,2),a0,f0,freq);
    
    % A(z,x) = az(z) .* ax(x), for discretizing  d/dz( A d/dz)
    az = @(z) 1./ez(z); ax = @(x) ex(x);
    
    % B(z,x) = bz(z) .* bx(x), for discretizing  d/dx( B d/dx)
    bz = @(z) ez(z); bx = @(x) 1./ex(x);
    
    C = @(z,x) ez(z) .* ex(x);
    
    iz = vec(1:nz); ix = vec(1:nx);
    [IZ,IX] = ndgrid(iz,ix);
    
    % circular shifts to compensate for the offsets created by spdiags
    Lz = @(xoffset) kron( spdiags(circshift(ax(ix+xoffset),xoffset),xoffset,nx,nx) , ...
                          spdiags( [circshift(az(iz-0.5),-1), ...
                                   -az(iz-0.5)-az(iz+0.5), ...
                                    circshift(az(iz+0.5),1)]/hz^2,-1:1,nz,nz) );
    
    Lx = @(zoffset) kron( spdiags( [circshift(bx(ix-0.5),-1), ...
                                   -bx(ix-0.5)-bx(ix+0.5), ...
                                    circshift(bx(ix+0.5),1)]/hx^2,-1:1,nx,nx),spdiags(circshift(bz(iz+zoffset),zoffset),zoffset,nz,nz) );
    
    d2z = b*Lz(0) + (1-b)/2 * (Lz(1) + Lz(-1));
    d2x = b*Lx(0) + (1-b)/2 * (Lx(1) + Lx(-1));
    
    Bk = sdiag(vec( C(IZ,IX) ));
    
    % Averaging matrices
    I0 = spdiags( 1/4*ones(N,4),[-nz -1 1 nz],N,N);
    I45 = spdiags( 1/4*ones(N,4),[-nz+1,-nz-1,1-nz,nz+1],N,N);    
    
    A = c*speye(N) + d*I0 + e*I45;
    
    L = d2z + d2x;
    H = L + A*sdiag(Bk*vec(fm));
    dH = A*sdiag(Bk*vec(dfm));
    if norm(ddfm)<1e-10, ddH = sdiag(zeros(size(ddfm)));         
    else ddH = A*sdiag(Bk*vec(ddfm)); end
end

function func = pml_func1d(nx,nb,a0,f0,freq)
% Inputs to func : dimensionless spatial coordinates
% func - 1d quadratic pml function
    dist_from_int = @(x) (x <= nb(1)) .* (nb(1)-x) + (x > nx-nb(2)) .* (x-(nx-nb(2)+1));
    sigma = @(x) 2*pi*a0*f0*( (dist_from_int(x)/max(nb)).^2 );
    func = @(x) 1-1i*sigma(x)/freq;
end



