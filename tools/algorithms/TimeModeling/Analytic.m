function [u] = Analytic(model,source,v)

%% Computes the 2D analytical solution for a constant velocity model and a given source at position xsrc,zsrc

% Ricker
t=1e-3*(0:model.dt:model.T);
[R,faxis] = fftrl(source,t,1);
%% 
nx=model.n(2);
nz=model.n(1);
h=model.d(1);

u=zeros(nx,nz,length(faxis));
for a=1:length(faxis)
    k = 2*pi*faxis(a)/v;
    for m=1:nx
        for n=1:nz
            tmp = k*sqrt((h*m - model.xsrc)^2 + (h*n - model.zsrc)^2);
            u(m,n,a)= -1i*pi*besselh(0,2,tmp)*R(a);
        end
    end
end

u=real(ifft(u,[],3));