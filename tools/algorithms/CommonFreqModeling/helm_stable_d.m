function d  = helm_stable_d( v, unit, scheme,freq )
%HELM_STABLE_D - Stable discretization gridpoint distance for the helmholtz
%equation
%
% Author:
%   Curt Da Silva, 2015
% 
%  Usage:
%    d = helm_stable_d( v, unit, scheme, freq )
%
%  Input:
%    v      - velocity model
%    unit   - one of 'm/s' (velocity) or 's2/m2' (slowness)
%    scheme - either PDEopts.HELM2D_CHEN9P, PDEopts.HELM3D_OPERTO27, PDEopts.HELM3D_STD7
%    freq   - maximum frequency 
% 
%  Output:
%    d      - 1 x 3 vector of computational domain spacings to use
   

    
switch scheme
  case PDEopts.HELM2D_CHEN9P
    npts_wavelength = 6;
  case PDEopts.HELM3D_OPERTO27
    npts_wavelength = 8;
  case PDEopts.HELM3D_CHEN27
    npts_wavelength = 8;
  case PDEopts.HELM3D_STD7
    npts_wavelength = 12;
  otherwise
    error(['Unknown scheme "', scheme, '"; Aborting...']);
end

switch unit
  case 's2/m2'
    vmin = min( v(:).^(-1/2));  
  case 's2/km2'
    vmin = min(1e3*(v(:).^(-1/2)));
  case 'm/s'
    vmin = min(v(:));
  otherwise
    error('Unknown unit');
end

d = [1 1 1]*vmin/(npts_wavelength*freq);


end

