function output = Ktomo(v0,input,flag,model)
% Scattering operator for traveltime tomography
%
% use:
%   output = Ktomo(v0,input,flag,model)
%
% input:
%   v0 - background velocity (scalar)
%   input - input vector
%   flag  - 1:forward, -1:adjoint
%   model - struct with parameters (see F.m)
%
% output:
%   output - result

%% 
N = prod(model.n);
[z,x]   = odn2grid(model.o,model.d,model.n);
[zz,xx] = ndgrid(z,x);

[zzs,xxs] = ndgrid(model.zsrc,model.xsrc);
[zzr,xxr] = ndgrid(model.zrec,model.xrec);

nrec  = length(model.xrec)*length(model.zrec);
nsrc  = length(model.xsrc)*length(model.zsrc);
nfreq = length(model.freq);


% define wavelet
omega = 2*pi*model.freq;
w = exp(1i*omega*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (omega).^2.*exp(-(omega/(2*pi*model.f0)).^2).*w;
end

% compute reference data and normalization factor
scale = zeros(nrec,nsrc);
Dref  = zeros(nrec,nsrc,nfreq);
for k = 1:nfreq
    Dref(:,:,k) = w(k)*prod(model.d)*Gconst(zzr(:),xxr(:),zzs(:),xxs(:),v0,omega(k));
    scale = scale + omega(k)^2*abs(Dref(:,:,k)).^2;
end

%% input
if flag > 0 
	output = zeros(nrec,nsrc);
else
	input  = reshape(input,[nrec,nsrc])./scale;
	output = zeros(N,1);
end

%% action of scattering operator
for k = 1:nfreq % loop over frequencies
	Gs = w(k)*prod(model.d)*Gconst(zz(:),xx(:),zzs(:),xxs(:),v0,omega(k));
    Gr = prod(model.d)*transpose(Gconst(zz(:),xx(:),zzr(:),xxr(:),v0,omega(k)));
	
	if flag > 0 % forward mode
		D_scat = omega(k)^2*Gr*diag(input)*Gs;
		output = output + 1i*omega(k)*conj(Dref(:,:,k)).*D_scat;
	else        % adjoint mode
		input_k = -1i*omega(k)*Dref(:,:,k).*input;
		output  = output + omega(k)^2*sum((Gr'*input_k).*conj(Gs),2);
	end
end

%% output
output = 1e-6*real(output);
if flag > 0
	output = vec(output./scale);
end



end


function g = Gconst(z,x,zs,xs,v0,omega)
% Greens function for constant velocity
%
% use: g = Gconst(z,x,zs,xs,v0,omega)
%
% input:
%   z,x   - grid
%   zs,xs - source locations
%   v0    - velocity (m/s)
%   omega - frequency
%
% output:
%   g - Greens function, each column is a source
%
    ns = length(zs);
    n  = length(z);
    g  = zeros(n,ns);
    for k = 1:length(zs)
        rr      = sqrt((x-xs(k)).^2 + (z-zs(k)).^2);  
        g(:,k)  = .25*1i*(besselj(0,omega*rr/v0) - 1i*bessely(0,omega*rr/v0)); g(x==xs(k)&z==zs(k),k) = -.5 + .25*1i;
    end
end