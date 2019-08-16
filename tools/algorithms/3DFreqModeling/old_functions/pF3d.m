function [D,J] = pF3d(m,Q,model)
% Frequency domain parallel 3D FD modeling operator
% use:
%   [D,J] = pF3d(m,Q,model)
% input:
%   m                 - vector with gridded squared slowness in [s^2/km^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition. size(Q,2) determines the number of
%                       sources
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc,ysrc} - vectors describing source array
%   model.{zrec,xrec,yrec} - vectors describing receiver array.
%
% output:
%   D  - Data cube (nrec x nsrc x nfreq) as vector. nsrc  = size(Q,2);
%                                                      nrec  = length(zrec)*length(xrec)
%                                                      nfreq = length(freq)
%   J  - Jacobian as SPOT operator  (still working on it)
% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


% pml
beta = 100;
pml = 2;

% physical grid
[z,x,y]  = odn2grid(model.o,model.d,model.n);

% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
[zt,xt,yt] = odn2grid(ot,dt,nt);
Nt = prod(nt);

% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec)*length(model.yrec);
nfreq  = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end

% if isfield(model,'wavelet')
% 	w = w.*model.wavelet;
% end

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(yt,model.yrec),opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(yt,model.ysrc),opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(3),model.nb(3)),opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));

% loop over frequencies
D = zeros([nrec,nsrc,nfreq]);
fid = fopen('F3d.log','w');
fprintf(fid,'freq, src, niter, res\n');
for k = 1:nfreq
    if model.np < 2
        [R,idx] = R_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml);
        a = 1./sqrt(sum(abs(R).^2,2));
        R = bsxfun(@times,R,a);
        for l = 1:nsrc
            opts.tol   = model.tol;
            opts.maxit = 5000;
            opts.w     = 1.5;
            [Uk,res]   = CGMN(R,idx,w(k)*a.*(Ps'*(Q(:,l))),zeros(Nt,1),opts);
            fprintf(fid,'%4d, %3d, %5d, %1.2e\n',k,l,length(res),res(end));
            D(:,l,k)   = Pr*Uk;
        end
        
    else
        mf      = 1e-6*Px*m;
        mf      = reshape(mf,nt);
        mdist   = distributed(mf);
        b       = ones(nt);
        b       = distributed(b);
        [R,idx] = pR_Helm3D(model.freq(k),mdist,b,ot,dt,nt,model.nb,beta,pml);
        spmd
            codist = getCodistributor(R);
            RLocal = getLocalPart(R);
            aLocal = 1./sqrt(sum(abs(RLocal).^2,2));
            for i = 1:27
                RLocal(:,i) = RLocal(:,i).*aLocal;
            end
            R           = codistributed.build(RLocal,codist,'noCommunication');
            Part        = codist.Partition;
            codist2     = codistributor1d(1, Part, [prod(nt),1]);
            a           = codistributed.build(aLocal,codist2,'noCommunication');
        end
        a1         = gather(a);
        
        spmd
            codist       = getCodistributor(R);
            Part         = codist.Partition;
            codist2      = codistributor1d(1, Part, [prod(nt),1]);
            x0           = codistributed.zeros([prod(nt),1],codist2,'noCommunication');
        end
        
        
        
        for l = 1:nsrc
            Ql         = Ps'*(Q(:,l));
            Ql         = w(k)*a1.*Ql;
            Ql         = distributed(Ql);
            spmd,
                Ql           = redistribute(Ql,codist2);
            end
            opts.tol   = model.tol;
            opts.maxit = 5000;
            opts.w     = 1.5;
            [Uk,res]   = pCARPCG(R,idx,Ql,x0,nt,opts);
            fprintf(fid,'%4d, %3d, %5d, %1.2e\n',k,l,length(res),res(end));
            Uk         = gather(Uk);
            D(:,l,k)   = Pr*Uk;
        end
    end
end
D = D(:);

if nargout > 1
   J = opppDF3d(m(:),Q,model);
end

