function [f,g] = misfit(m,Q,Is,eps,D,model,params)
% misfit for 3D FWI using CGMN.
%
% use:
%   [f,g] = misfit(m,Q,Is,eps,D,model,params)
%
% input:
%   m   - model
%   Q   - source
%   Is  - source indices to be used
%   eps - tolerance for forward and adjoint solves 
%   D   - data
%   model - struct with model parameters, see F.m
%   params.srcest - estimate source weight {0,1}
%   params.hmin   - min. offset to use
%   params.maxit  - max. iterations for CGMN solvess
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
D = reshape(D,[nrec nsrc nfreq]);
f = 0;
g = zeros(N,1);

opts.tol   = eps;
opts.maxit = params.maxit;
opts.w     = 1.5;
fid = fopen('misfit.log','w');
fprintf(fid,'freq, src, niter, res\n');

for k = 1:nfreq
    [R,idx] = R_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml);
    a = 1./sqrt(sum(abs(R).^2,2));
    G = G_Helm3D(model.freq(k),1e-6*(Px*m),ones(Nt,1),ot,dt,nt,model.nb,beta,pml); 
    for l = Is
        % forward equation
        [U0k,res]   = CGMN(bsxfun(@times,R,a),idx,w(k)*a.*(Ps'*(Q(:,l))),zeros(Nt,1),opts);
        fprintf(fid,'%4d, %3d, %5d, %1.2e\n',k,l,length(res),res(end));
        
        % predicted data
        Dpred = Pr*U0k;
        
        % offset weight
        Wh = offset_weight(params.hmin,l,model);
        
        % source weight
        if params.srcest
            wkl = ((Wh.*Dpred)'*(Wh.*D(:,l,k)))/norm(Wh.*Dpred).^2;
        else
            wkl = 1;
        end
        
        % local misfit
        fkl      = .5*norm(wkl*Wh.*Dpred - Wh.*D(:,l,k))^2;
        
        % total misfit
        f = f + fkl;
        
        if nargout > 1
            % adjoint equation
            [R,idx]   = Rtransp(R,idx);
            a = 1./sqrt(sum(abs(R).^2,2));
            
            [V0k,res] = CGMN(conj(bsxfun(@times,R,a)),idx,a.*(Pr'*(wkl'*Wh.*(wkl*Wh.*Dpred - Wh.*D(:,l,k)))),zeros(Nt,1),opts);
            fprintf(fid,'%4d, %3d, %5d, %1.2e\n',k,l,length(res),res(end));

            % gradient
            g = g - 1e-6*Px'*(real(conj(U0k).*(G'*V0k)));
            
            % transpose back
            [R,idx] = Rtransp(R,idx);
            a = 1./sqrt(sum(abs(R).^2,2));
        end
        
    end
end
fclose(fid);

% scale misfit and  gradient
f = f/length(Is);
g = g/length(Is);

end

function Wh = offset_weight(hmin,l,model)
    [iz,ix,iy]=ind2sub([length(model.zsrc) length(model.xsrc) length(model.ysrc)],l);
    [hz,hx,hy]=ndgrid(model.zrec-model.zsrc(iz),model.xrec-model.xsrc(ix),model.yrec-model.ysrc(iy));
    Wh = sqrt((length(model.yrec)>1)*hy.^2 + (length(model.xrec)>1)*hx.^2+(length(model.zrec)>1)*hz.^2)>hmin;
    Wh=Wh(:);
end