function w = fwi_wavelet(f,t0,f0,type)
% FWI_WAVELET
%   Generates a wavelet signature for FWI simulations
%   
% Usage:
%   w = fwi_wavelet(f,t0,f0,type);
%  
% Input:
%   f    - frequency grid (vector)
%   t0   - firing time
%   f0   - if specified, peak frequency of the Ricker wavelet, 
%          otherwise a wavelet with a flat spectrum is generated
%   type - one of
%          PDEopts.WAVELET_RICKER  - Ricker wavelet

% define wavelet
if exist('type','var')==0
    type = PDEopts.WAVELET_RICKER;
end
   
switch type
    case PDEopts.WAVELET_RICKER
      
      w = exp(1i*2*pi*f*t0);
      if exist('f0','var') && ~isempty(f0) && f0
          % Ricker wavelet with peak-frequency model.f0
          w = f.^2 .* exp(-(f./f0).^2) .* w;
      end
    
end