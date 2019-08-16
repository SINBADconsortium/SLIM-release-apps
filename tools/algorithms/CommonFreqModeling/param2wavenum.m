function [ v, dv, ddv ] = param2wavenum( v, freq, unit )
%PARAM2WAVENUM Converts the input parameter to wavenumber for use in the
%discrete_helmholtz function. 
%
%  Curt Da Silva, 2015
%
%  Usage:
%    [v,dv,ddv] = param2wavenum(v, freq, unit);
%
%  Input:
%    v     - velocity/slowness^2 parameter
%    freq  - frequency
%    unit  - one of 'm/s','s2/m2', or 's2/km2'
% 
%  Output:
%    v     - wavenumber parameter (rad^2 * s^2/m^2)
%    dv    - first derivative of v
%    ddv   - second derivative of v
%
c = (2*pi*freq)^2;
switch unit
  case 'm/s'
    if nargout >=2
        dv = -2*c*(v.^(-3));
    end
    if nargout >=3
        ddv = 6 * c*(v.^(-4));
    end
    v = c*(v.^(-2));
    
  case 's2/m2';
    v = c*v;
    if nargout >=2
        dv = c*ones(size(v)); 
    end
    if nargout >=3
        ddv = zeros(size(v)); 
    end
  case 's2/km2';
    if nargout >= 2
        dv = 1e-6*c*ones(size(v));
    end
    if nargout >=3
        ddv = zeros(size(v));
    end
    v = 1e-6*c*v;
    otherwise
        error(['Unknown unit - ' unit]);
end

