function output = DF_4D(m,Q,input,flag,model,dogather)
% Frequency domain modeling in the Born approximation. This is the
% Jacobian of F_4D(m,Q,model). 
%
% use: 
%   output = DF_4D(m,Q,input,flag,model,{gather})
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source and has to be
%                       distributed over the last dimension.
%   input             - flag= 1: vector with gridded slowness perturbation
%                       flag=-1: vectorized data cube of size nrec xnrec x nfreq
%   flag              -  1: forward mode
%                       -1: adjoint mode
%   model.{o,d,n}     - regular grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc
%   model.nb          - number of extra points for absorbing boundary on each side
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array
%   gather            - indicate whether to gather the output (optional, default = 0)
%
%%
% *Example*
% Below example defines a grid [0,1000]^2 with 10 m gridspacing. The
% velocity is 2 km/s. The sources and receivers coincide at x = [0:10:1000]
% and z = 10m. We do the dottest with a randomly generated perturbation.
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
% dm = randn(size(m));
% dD = DF_4D(m,Q,dm,1,model);
% dDt = 0*dD + randn(size(dD));
% dmt = DF_4D(m,Q,dDt,-1,model);
% real(dD'*dDt)
% dm'*dmt
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

if nargin < 6
    dogather = 0;
end

% handle 3rd dimension consistently
model.n(3) = 1; model.o(3) = 0;model.nb(3) = 0; model.d(3) = 1;

% physical grid
[z,x]  = odn2grid(model.o,model.d,model.n);

% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
Nt = prod(nt);
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
Pe = opKron(opExtension(model.n(2),model.nb(2),0),opExtension(model.n(1),model.nb(1),0));
% model parameter
nu  = 1e-3*sqrt(Px*m);
dnu = opDiag(.5*1e-6./nu);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);
w    = distributed(w);
if flag==1
    % solve Helmholtz for each frequency in parallel
    spmd
        codistr   = codistributor1d(2,codistributor1d.unsetPartition,[nsrc*nrec,length(freq)]);
        freqloc   = getLocalPart(freq);
        wloc      = getLocalPart(w);
        nfreqloc  = length(freqloc);
        outputloc = zeros(nsrc*nrec,nfreqloc);
        if size(Q,3) == 1
            for k = 1: nfreqloc
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Q));
                Sk        = -(2*pi*freqloc(k))*(dnu*(dHk*(U0k.*repmat(Px*input,1,nsrc))));
                U1k       = Hk\Sk;
                outputloc(:,k) = vec(Pr*U1k);
            end
        else
            Qloc = getLocalPart(Q);
            for k = 1: nfreqloc
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Qloc(:,:,k)));
                Sk        = -(2*pi*freqloc(k))*(dnu*(dHk*(U0k.*repmat(Px*input,1,nsrc))));
                U1k       = Hk\Sk;
                outputloc(:,k) = vec(Pr*U1k);
            end
        end
        output = codistributed.build(outputloc,codistr,'noCommunication'); 
    end
    output = vec(output);
    if dogather
        output = gather(output);
    end
else
    if ~isdistributed(input)
        input = distributed(input);
    end
    input = invvec(input,[nrec*nsrc,nfreq]);
    spmd
        freqloc   = getLocalPart(freq);
        wloc      = getLocalPart(w);
        nfreqloc  = length(freqloc);
        outputloc = zeros(prod(model.n),1);
        inputloc  = getLocalPart(input);
        if size(Q,3)==1
            for k = 1:nfreqloc
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Q));
                Sk        = -Pr'*reshape(inputloc(:,k),[nrec nsrc]);
                V0k       = Hk'\Sk;
                r         = (2*pi*freqloc(k))*real(sum(conj(U0k).*(dHk'*(dnu'*V0k)),2)); 
                outputloc = outputloc + Pe'*r;
            end
        else
            Qloc = getLocalPart(Q);
            for k = 1:nfreqloc
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Qloc(:,:,k)));
                Sk        = -Pr'*reshape(inputloc(:,k),[nrec nsrc]);
                V0k       = Hk'\Sk;
                r         = (2*pi*freqloc(k))*real(sum(conj(U0k).*(dHk'*(dnu'*V0k)),2)); 
                outputloc = outputloc + Pe'*r;
            end
        end
        output = pSPOT.utils.global_sum(outputloc);
    end
    output = output{1};
end

