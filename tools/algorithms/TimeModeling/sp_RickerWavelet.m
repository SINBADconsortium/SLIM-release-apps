function s=sp_RickerWavelet(f0,t0,dt,T)

%function S=SP_RICKERWAVELET(F0,T0,DT,T)
%
% Ricker wavelet of peak frequency f0, scaled so as to have a maximum
% amplitude of 1, and shifted of t0 in time.
% The scale factor, alpha, is given by the expression:
%
%   alpha = sqrt(3*sigma)*pi^.25 / 2 , where sigma = 1/(sqrt(2)*pi*f0)
%
% INPUTS
%
%   f0 [Hz]: peak frequency
%   t0 [s]: shift in time
%   dt [s]: time step
%   T [s]: acquisition duration
% Author : Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015

t = (0:dt:T);
r = (pi*f0*(t-t0));
q = (1-2.*r.^2).*exp(-r.^2);
s=q';

end