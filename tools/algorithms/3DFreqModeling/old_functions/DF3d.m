function output = DF3d(m,Q,input,flag,model)
% Frequency domain 3D FD linearized modeling operator
%
% use: 
%   output = DF3d(m,Q,input,flag,model)
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition. size(Q,2) determines the number of
%                       sources
%   input             - model perturbation or data perturbation
%   flag              - 1: foward mode, -1: adjoint mode
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc,ysrc} - vectors describing source array
%   model.{zrec,xrec,yrec} - vectors describing receiver array.
%   model.tol         - tolerance for CGMN
%
% output:
%   output - data perturbation or model perturbation
%
% Author: Tristan van Leeuwen
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
N = prod(model.n);

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

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(yt,model.yrec),opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(yt,model.ysrc),opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(3),model.nb(3)),opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));


% loop over frequencies
if flag == 1
    output = zeros(nrec,nsrc,nfreq);
    for k = 1:nfreq
        [R,idx] = R_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml);
        a = 1./sqrt(sum(abs(R).^2,2));
        R = bsxfun(@times,R,a);
        G = G_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml);
        for l = 1:nsrc
            opts.tol   = model.tol;
            opts.maxit = 5000;
            opts.w     = 1.5;
            U0k        = CGMN(R,idx,w(k)*a.*(Ps'*(Q(:,l))),zeros(Nt,1),opts);
            U1k        = CGMN(R,idx,-1e-6*a.*(G*(U0k.*(Px*input))),zeros(Nt,1),opts);
            output(:,l,k) = Pr*U1k;
        end
    end
    output = output(:);
else
    input  = reshape(input,nrec,nsrc,nfreq);  
    output = zeros(N,1);
    for k = 1:nfreq
        [R,idx] = R_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml);
        a = 1./sqrt(sum(abs(R).^2,2));
        G       = G_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml);
        for l = 1:nsrc
            opts.tol   = model.tol;
            opts.maxit = 5000;
            opts.w     = 1.5;
            
            % forward
            U0k        = CGMN(bsxfun(@times,R,a),idx,w(k)*a.*(Ps'*(Q(:,l))),zeros(Nt,1),opts);
            
            % transpose
            [R,idx]    = Rtransp(R,idx);
            a = 1./sqrt(sum(abs(R).^2,2));
            
            % backward
            V0k        = CGMN(bsxfun(@times,conj(R),a),idx,a.*(Pr'*input(:,l,k)),zeros(Nt,1),opts);
            
            % transpose back
            [R,idx] = Rtransp(R,idx);
            a = 1./sqrt(sum(abs(R).^2,2));
            
            % output
            output   = output - 1e-6*Px'*(real(conj(U0k).*(G'*V0k)));
        end
    end
    output = output(:);
end



