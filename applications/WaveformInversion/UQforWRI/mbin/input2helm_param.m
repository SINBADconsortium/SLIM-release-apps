function [f,df,ddf] = input2helm_param(m,freq)
% INPUT2HELM_PARAM - Converts the s^2/km^2 input to  
% s/(m) * radians. The output of this function is the input to Helm2D
%
% Usage:
%   [f,df,ddf] = input2helm_param(m,freq);
%  
% Input:
%   m    - model parameter
%   freq - frequency
%  
% Output:
%   f    - model parameters in s/(km) *radians
%   df   - gradient of f
%   ddf  - hessian of f
f = 2*pi*freq* 1e-3 *(m).^(1/2);
df = 2*pi*freq* 1e-3 * ( 0.5  * (m.^(-1/2)) );
ddf = 2*pi*freq*1e-3 * (- 0.25 * (m.^(-3/2) ) );
end


