function D = G(m,Q,model)
% 2D Greens function for Helmholtz equation with linear velocity profile 
% v = v0 + alpha*z. NOT STABLE FOR SMALL (<.2) NON-ZERO VALUES OF ALPHA!
%
% The output is consistent with that of F.m. 
% Usefull for testing the accuracy of F.m, subtracting direct wave or
% finding linear starting model.
%
% use: 
%   [D] = G(m,Q,model)
% input:
%   m                 - [v0;alpha] defining linear velocity profile v0+alpha*z
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source and has to be
%                       distributed over the last dimension.
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%
% output:
%   D  - Data cube (nsrc x nrec x nfreq) as vector. nsrc  = size(Q,2);
%                                                   nrec  = length(zrec)*length(xrec) 
%                                                   nfreq = length(freq)

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

v0 = m(1);
alpha = m(2);

[zz,xx]   = ndgrid(model.zrec,model.xrec);
[zzs,xxs] = ndgrid(model.zsrc,model.xsrc);

nrec  = length(model.xrec)*length(model.zrec);
nsrc  = length(model.xsrc)*length(model.zsrc);
nfreq = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end

D = zeros(nsrc*nrec,nfreq);

for k = 1:nfreq
    % data
    Dk = zeros(nrec,nsrc);
    for l = 1:nsrc
        if alpha==0
            G = Gconst(zz,xx,zzs(l),xxs(l),v0,model.freq(k));
        else
            G = Glin(zz,xx,zzs(l),xxs(l),v0,alpha,model.freq(k));
        end
        Dk(:,l) = w(k)*vec(G);
    end
    % multiply with source function.
    if size(Q,3)==1
        D(:,k) = model.d(1)*model.d(2)*vec(Dk*Q);
    else
        D(:,k) = model.d(1)*model.d(2)*vec(Dk*Q(:,:,k));
    end
end

D = vec(D);

end

function g = Gconst(z,x,zs,xs,v0,f)

rr = sqrt((x-xs).^2 + (z-zs).^2);
w  = 2*pi*f;
g  = .25*1i*(besselj(0,w*rr/v0) - 1i*bessely(0,w*rr/v0)); g(x==xs&z==zs) = -.5 + .25*1i;

end

function g = Glin(z,x,zs,xs,v0,alpha,f)
% see Kuvshinov & Mulder: The exact solution of the time-harmonic wave equation for a linear
% velocity profile. GJI (2006) 167, 659–662.

nu = 1i*sqrt((2*pi*f/alpha).^2 - 0.25);

u = 1+((x-xs).^2+(z-zs).^2)./(2*(zs+v0/alpha).*(z+v0/alpha));

g = -(1/(2*pi))*legendreQ(nu-.5,u);

end



