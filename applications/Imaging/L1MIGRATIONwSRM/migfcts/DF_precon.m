function [output, precon] = DF_precon(m,Q,input,flag,model)
% Frequency domain modeling by preconditioned Born approximation. Weighted by
% the inverse of incident wavefield square.
%
% use: 
%   [output, precon] = DF_precon(m,Q,input,flag,model)
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
%	model.winsize     - length of smooth window, a vector consistent with o/d/n
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
% model.freq = [5]; 
% model.f0 = 0;
% model.t0 = 0;
% model.zsrc = 10; 
% model.xsrc = 100;
% model.zrec = 10; 
% model.xrec = 0:10:1000;
% model.winsize = [500,500,1];
% m = .25*ones(prod(model.n),1);
% Q = speye(length(model.xsrc));
% dm = randn(size(m));
% dD = DF_precon(m,Q,dm,1,model);
% dDt = 0*dD + randn(size(dD));
% dmt = DF_precon(m,Q,dDt,-1,model);
% real(dD'*dDt)
% dm'*dmt
%%

% Author: Ning Tu (Adapted from Tristan's DF function)
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: April, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

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

% model parameter
nu  = 1e-3*sqrt(Px*m);
dnu = opDiag(.5*1e-6./nu);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);
w    = distributed(w);

% check source matrix input
if (size(Q,3)==1)&&(isdistributed(Q))
    Q = gather(Q);
end

% smoothing parameters for preconditioning
kSmooth = ceil(model.winsize./model.d);
kSmooth = ceil((kSmooth-1)/2);
opDoSmooth = opKron(opDirac(nsrc),opKron(opSmooth(nt(2),kSmooth(2)),opSmooth(nt(1),kSmooth(1))));

if flag==1
    % solve Helmholtz for each frequency in parallel
    spmd
        codistr   = codistributor1d(2,codistributor1d.unsetPartition,[nsrc*nrec,length(freq)]);
        freqloc   = getLocalPart(freq);
        wloc      = getLocalPart(w);
        nfreqloc  = length(freqloc);
        outputloc = zeros(nsrc*nrec,nfreqloc);
        preconloc = zeros(nrec*nsrc,nfreqloc);
        if size(Q,3) == 1
            for k = 1: nfreqloc
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Q));
                EU0k      = opDoSmooth*real(vec(conj(U0k).*(U0k)));
                if not(all(EU0k))
                    error(['Fatal: devided by 0. Increase the smooth diameter.']);
                end
                EU0k      = invvec(EU0k,[Nt,nsrc]);
                Sk        = -(2*pi*freqloc(k))*(dnu*(dHk*(U0k.*repmat(Px*input,1,nsrc)./EU0k)));
                U1k       = Hk\Sk;
                preconloc(:,k) = vec(Pr*EU0k);
                outputloc(:,k) = vec(Pr*U1k);
            end
        else
            Qloc = getLocalPart(Q);
            for k = 1: nfreqloc
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Qloc(:,:,k)));
                EU0k      = opDoSmooth*real(vec(conj(U0k).*(U0k)));
                if not(all(EU0k))
                    error(['Fatal: devided by 0. Increase the smooth diameter.']);
                end
                EU0k      = invvec(EU0k,[Nt,nsrc]);
                Sk        = -(2*pi*freqloc(k))*(dnu*(dHk*(U0k.*repmat(Px*input,1,nsrc)./EU0k)));
                U1k       = Hk\Sk;
                preconloc(:,k) = vec(Pr*EU0k);
                outputloc(:,k) = vec(Pr*U1k);
            end
        end
        output = codistributed.build(outputloc,codistr,'noCommunication'); 
        precon = codistributed.build(preconloc,codistr,'noCommunication'); 
    end
    output = vec(output);
else
    spmd
    	codistr   = codistributor1d(2,codistributor1d.unsetPartition,[nsrc*nrec,length(freq)]);
        freqloc   = getLocalPart(freq);
        wloc      = getLocalPart(w);
        nfreqloc  = length(freqloc);
        outputloc = zeros(prod(model.n),1);
        inputloc  = getLocalPart(input);
        preconloc = zeros(nrec*nsrc,nfreqloc);
        if size(Q,3)==1
            for k = 1:nfreqloc
                inputloc  = reshape(inputloc,[nsrc*nrec,nfreqloc]);
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Q));
                Sk        = -Pr'*reshape(inputloc(:,k),[nrec nsrc]);
                V0k       = Hk'\Sk;
                EU0k      = opDoSmooth*real(vec(conj(U0k).*(U0k)));
                if not(all(EU0k))
                    error(['Fatal: devided by 0. Increase the smooth diameter.']);
                end
                EU0k      = invvec(EU0k,[Nt,nsrc]);
                r         = (2*pi*freqloc(k))*real(sum(conj(U0k).*(dHk'*(dnu'*V0k))./EU0k,2));
                preconloc(:,k) = vec(Pr*EU0k);
                outputloc = outputloc + Px'*r;
            end
        else
            Qloc = getLocalPart(Q);
            for k = 1:nfreqloc
                inputloc  = reshape(inputloc,[nsrc*nrec,nfreqloc]);
                [Hk, dHk] = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
                U0k       = Hk\(wloc(k)*(Ps'*Qloc(:,:,k)));
                Sk        = -Pr'*reshape(inputloc(:,k),[nrec nsrc]);
                V0k       = Hk'\Sk;
                EU0k      = opDoSmooth*real(vec(conj(U0k).*(U0k)));
                if not(all(EU0k))
                    error(['Fatal: devided by 0. Increase the smooth diameter.']);
                end
                EU0k      = invvec(EU0k,[Nt,nsrc]);
                r         = (2*pi*freqloc(k))*real(sum(conj(U0k).*(dHk'*(dnu'*V0k))./EU0k,2));
                preconloc(:,k) = vec(Pr*EU0k);
                outputloc = outputloc + Px'*r;
            end
        end
        output = pSPOT.utils.global_sum(outputloc);
        precon = codistributed.build(preconloc,codistr,'noCommunication'); 
    end
    output = output{1};
end


