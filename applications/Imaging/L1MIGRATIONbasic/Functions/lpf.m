%function [f] = filter2d(A,alfa)
% Returns a 2D filtered version of A. 
function [f] = lpf(X,alfa)
X=extender2(X);

[M N] = size(X);

kk=[-M/2+1:M/2]';
jj=[-N/2+1:N/2];


b=alfa*N/2;
a=alfa*M/2;

% size(kk)
% size(jj)

fil=fftshift(exp(-(repmat(kk.^2,1,N)/a^2)).*exp(-(repmat(jj.^2,M,1)/b^2)));


% mesh(fftshift(fil)); colorbar
% pause(1);
f=real(ifft2(fft2(X).*fil));
f=extender2(f,'transp');

