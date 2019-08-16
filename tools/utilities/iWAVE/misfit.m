function [f, g] = misfit(m,model,D,options)
%% Misfit function for FWI
% Usage:
% [f, g] = misfit(m,model,D,options)
%
% Input:
%   m                 - vector with gridded parameters [bulk [10^9 kg/(m*s^2)]; buoyancy [cm^3/gr]]
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%   options           - options for iWAVE
%   options.fwdpara   - the forward operator parameter file
%   options.adjpara   - the adjoint operator parameter file
%   options.linpara   - the linear operator parameter file
%   options.delete      - whether delete or not delete the file that iWave creates
%
% output:
%   f  - misfit
%   g  - gradient
%
%
% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: August, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.





[D1,J]     = Fd2DT(m,model,options);
f          = 0.5*norm(D1(:)-D(:))^2;
if nargout > 1
	g      = J' * (D1(:)-D(:));
	if options.bulkonly == 1
		g = [g;zeros(prod(model.n),1)];
	end
end

