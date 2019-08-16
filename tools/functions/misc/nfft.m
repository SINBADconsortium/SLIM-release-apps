% NFFT.M
% 
%  Computes the nonequispaced fft and its adjoint.
% 
%  f=nfft(f_hat,x,options,sigma,m)
% 
%  f_hat   Fourier coefficients
%  x       nodes \subset [-1/2,1/2)
%  N       polynomial degree
% 
%  ft_opt  .method   'gaussian'     convolution based, using a gaussian
%                                   window function 
%                    'f_gaussian'   avoids direct calls of exp()
%                    'taylor'       taylor expansion on the nearest
%                                   neighbor on a fine grid
%          .sigma         oversampling factor
%          .m             cut-off or Taylor-order
% 
%  tflag   'notransp' evaluate trigonometric polynomial
%          'transp'   adjoint transform
% 
%  Author  Stefan Kunis
% 
% ---------------------------------------------------------------
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.
%
function b=nfft(a,x,N,ft_opt,tflag)
%
% Computes the nonequispaced fft and its adjoint.
%
% f=nfft(f_hat,x,options,sigma,m)
%
% f_hat   Fourier coefficients
% x       nodes \subset [-1/2,1/2)
% N       polynomial degree
%
% ft_opt  .method   'gaussian'     convolution based, using a gaussian
%                                  window function 
%                   'f_gaussian'   avoids direct calls of exp()
%                   'taylor'       taylor expansion on the nearest
%                                  neighbor on a fine grid
%         .sigma         oversampling factor
%         .m             cut-off or Taylor-order
%
% tflag   'notransp' evaluate trigonometric polynomial
%         'transp'   adjoint transform
%
% Author  Stefan Kunis

options=(ft_opt.method);
sigma=ft_opt.sigma;
m=ft_opt.m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonequispaced fft
%
%         N/2-1  
% f(j) =  sum    f_hat(k+N/2+1)*exp(-2*pi*i*k*x(j)), 1 <= j <= M.
%        k=-N/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(tflag,'notransp')

f_hat=a;
M=length(x);
n=sigma*2^ceil(log2(N));
sigma=n/N;

f=zeros(size(x));

switch options
  case 'fast_gaussian'
    freq=(-(N/2):(N/2-1))';
    b=2*sigma*m / ((2*sigma-1)*pi);
    inv_phi_hat=exp(b*(pi*freq/n).^2);
    g_hat=[zeros((n-N)/2,1);f_hat.*(inv_phi_hat);zeros((n-N)/2,1)];
    g=fft(fftshift(g_hat));
    exp_b=exp(-(0:2*m+1).^2 /b);
    for j=1:M
      c_j=n*x(j);
      u_j=ceil(c_j-m);
      o_j=floor(c_j+m);
      supp_j=mod(u_j:o_j,n)+1;
      psi_j=zeros(size(u_j:o_j));
      psi_j(1)=(pi*b)^(-1/2) * exp(-(n*x(j)-(u_j)).^2 /b);
      exp_x=exp(2*(n*x(j)-(u_j)) /b);
      for l=2:(o_j-u_j+1)
        psi_j(l)=psi_j(l-1) * exp_x;
      end;
      psi_j = psi_j.* exp_b(1:(o_j-u_j+1));
      f(j)=psi_j*g(supp_j);
    end;
  case 'gaussian'
    freq=(-(N/2):(N/2-1))';
    b=2*sigma*m / ((2*sigma-1)*pi);
    inv_phi_hat=exp(b*(pi*freq/n).^2);
    g_hat=[zeros((n-N)/2,1);f_hat.*(inv_phi_hat);zeros((n-N)/2,1)];
    g=fft(fftshift(g_hat));
    for j=1:M
      c_j=n*x(j);
      u_j=ceil(c_j-m);
      o_j=floor(c_j+m);
      supp_j=mod(u_j:o_j,n)+1;
      psi_j=(pi*b)^(-1/2) * exp(-(n*x(j)-(u_j:o_j)).^2 /b);
      f(j)=psi_j*g(supp_j);
    end;
  case 'taylor'
    freq=-2*pi*i*(-(N/2):(N/2-1))';
    ix=floor(n*(x+0.5)+1);
    dx=x-((ix-1)/n-0.5);
    for l=0:m
      g_hat=[zeros((n-N)/2,1);f_hat.*(freq.^l);zeros((n-N)/2,1)];
      g=fftshift(fft(fftshift(g_hat)));
      f=f+g(ix).*(dx.^l)/prod(1:l);
    end;
  otherwise
    disp('unknown option')
end;

b=f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjoint nonequispaced fft
%
%                    M  
% f_hat(k+N/2+1) =  sum    f(j)*exp(2*pi*i*k*x(j)), -N/2 <= k <= N/2-1.
%                   j=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(tflag,'transp')

f=a;
M=length(x);
n=sigma*2^ceil(log2(N));
sigma=n/N;

switch options
  case 'fast_gaussian'
    freq=(-(N/2):(N/2-1))';
    b=2*sigma*m / ((2*sigma-1)*pi);
    exp_b=exp(-(0:2*m+1).^2 /b);
    g=zeros(n,1);
    for j=1:M
      c_j=n*x(j);
      u_j=ceil(c_j-m);
      o_j=floor(c_j+m);
      supp_j=mod(u_j:o_j,n)+1;
      psi_j=zeros(size(u_j:o_j));
      psi_j(1)=(pi*b)^(-1/2) * exp(-(n*x(j)-(u_j)).^2 /b);
      exp_x=exp(2*(n*x(j)-(u_j)) /b);
      for l=2:(o_j-u_j+1)
        psi_j(l)=psi_j(l-1) * exp_x;
      end;
      psi_j = psi_j.* exp_b(1:(o_j-u_j+1));
      g(supp_j)=g(supp_j) + f(j)*psi_j';
    end;
    inv_phi_hat=exp(b*(pi*freq/n).^2);
    g_hat=n*fftshift(ifft(g));
    f_hat=g_hat((n/2-N/2+1):(n/2+N/2)).*(inv_phi_hat);
  case 'gaussian'
    freq=(-(N/2):(N/2-1))';
    b=2*sigma*m / ((2*sigma-1)*pi);
    g=zeros(n,1);
    for j=1:M
      c_j=n*x(j);
      u_j=ceil(c_j-m);
      o_j=floor(c_j+m);
      supp_j=mod(u_j:o_j,n)+1;
      psi_j=(pi*b)^(-1/2) * exp(-(n*x(j)-(u_j:o_j)).^2 /b);
      g(supp_j)=g(supp_j) + f(j)*psi_j';
    end;
    inv_phi_hat=exp(b*(pi*freq/n).^2);
    g_hat=n*fftshift(ifft(g));
    f_hat=g_hat((n/2-N/2+1):(n/2+N/2)).*(inv_phi_hat);
  case 'taylor'
    freq=-2*pi*i*(-(N/2):(N/2-1))';
    ix=floor(n*(x+0.5)+1);
    dx=x-((ix-1)/n-0.5);
    f_hat=zeros(N,1);
    for l=0:m
      g=[accumarray(ix,dx.^l .*f);zeros(n-max(ix),1)];
      g_hat=n*fftshift(ifft(fftshift(g)));
      f_hat=f_hat+(-freq).^l/prod(1:l).*g_hat((n/2-N/2+1):(n/2+N/2));
    end;
  otherwise
    disp('unknown option')
end;

b=f_hat;
end;
