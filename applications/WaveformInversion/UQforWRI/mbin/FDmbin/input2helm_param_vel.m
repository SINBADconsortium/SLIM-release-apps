function [f,df,ddf] = input2helm_param_vel(m,freq)
% INPUT2HELM_PARAM - Converts the s^2/m^2 input to  
% s/(km) * radians. The output of this function is the input to Helm2D
%
% Usage:
%   [f,df,ddf] = input2helm_param(m,freq);
%  
% Input:
%   m    - model parameter (velcity m/s)
%   freq - frequency
%  
% Output:
%   f    - model parameters in s/(km) *radians
%   df   - gradient of f
%   ddf  - hessian of f
f = 2*pi*freq./(m);
df = 2*pi*freq* (- 1./ m.^2);
ddf = 2*pi*freq * (2 ./ m.^3 );
end
