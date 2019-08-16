function [D,J] = F_4D(m,Q,model,dogather)
% Frequency domain FD modeling operator
%
% use: 
%   [D,J] = F_4D(m,Q,model,{gather})
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source and has to be
%                       distributed over the last dimension.
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%   gather            - 1: gather output, default = 0;
%
% output:
%   D  - Data cube (nrec x nsrc x nfreq) as (distributed) vector. nsrc  = size(Q,2);
%                                                                 nrec  = length(zrec)*length(xrec) 
%                                                                 nfreq = length(freq)
%   J  - Jacobian as pSPOT operator
%%
% *Example*
% Below example defines a grid [0,1000]^2 with 10 m gridspacing. The
% velocity is 2 km/s. The sources and receivers coincide at x = [0:10:1000]
% and z = 10m.
%
%%
% model.o = [0 0 0];
% model.d = [10 10 1];
% model.n = [101 101 1];        
% model.nb = [10 10 0];
% model.freq = [5 10 15 20]; 
% model.f0 = 0;
% model.zsrc = 10; 
% model.xsrc = 0:10:1000;
% model.zrec = 10; 
% model.xrec = 0:10:1000;
% m = .25*ones(prod(model.n),1);
% Q = speye(length(model.xsrc));
% D = F(m,Q,model);
% 
%%

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

if nargin < 4
    dogather = 0;
end

% handle 3rd dimension consistently
model.n(3)  = 1; model.o(3)  = 0; model.d(3)  = 1;
model.nb(3) = 0;

% physical grid
[z,x]  = odn2grid(model.o,model.d,model.n);

% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));

% model parameter: slowness [s/m] on computational grid.
nu = 1e-3*Px*sqrt(m);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);
w    = distributed(w);

spmd
    codistr  = codistributor1d(2,[],[nsrc*nrec,length(freq)]);
    freqloc  = getLocalPart(freq);
    wloc     = getLocalPart(w);
    nfreqloc = length(freqloc);
    Dloc     = zeros(nrec*nsrc,nfreqloc);
    if size(Q,3)==1
        for k = 1:nfreqloc
            Hk  = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
            Uk  = Hk\(wloc(k)*(Ps'*Q)); 
            Dloc(:,k) = vec(Pr*Uk); 
        end
    else
        Qloc = getLocalPart(Q);
        for k = 1:nfreqloc
            Hk  = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
            Uk  = Hk\(wloc(k)*(Ps'*Qloc(:,:,k))); 
            Dloc(:,k) = vec(Pr*Uk); 
        end
    end
    D = codistributed.build(Dloc,codistr,'noCommunication'); 
end

% vectorize output, gather if needed
D = vec(D);

if dogather 
    D = gather(D); 
end

% construct pSPOT operator
J = oppDF_4D(m,Q,model,dogather);
