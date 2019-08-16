


function [f g] = funSTst(r, params)


N = length(r);

para.n = N;
para.r = r;
para.k = params.nu;
para.scaleIn = 0.5;

a = ComputeScale(para);

%a = 0.01; %10*norm(r, inf);  % this is craaaazy
r = r./a;
rsq = conj(r).*r;  % rescaling 

% self tuning student's t function 

fh   = @(k)g_df(rsq,k,a,1);
%K = scnewton(fh,.5,1e-4);
 
K = scnewton(fh,params.nu,1e-4);
params.nu = K;

%if flag==1 % true
ftrue = @(k)(-N*log(gamma((k+1)/2)/(gamma(k/2)*sqrt(pi*k))) + ((k+1)/2)*sum(log(1+rsq/k)));
%h = 1e-3;
%f = fh(K);
f = ftrue(K);
g = ((K+1)/2).*r./(K+rsq);
%g = (ftrue(K+h) - ftrue(K-h))/(2*h);
%h = (ftrue(K+h) - 2*ftrue(K) + ftrue(K-h))/(h.^2);
    %h = N/(2*K*(K+1)^2) + (2*K-1)/(4*K.^2)*sum(R./(K+R)) + (K+1)/(2*K)*sum(R./(K+R).^2);
%else % approximation
%    f = (K*N/2)*log(K/(K+1)) + ((K+1)/2)*sum(log(1+R/K)) + (N/2)*(1+log(2*pi));
%    g = (N/2)*(log(K/(K+1)) + 1/(K+1)) + (1/2)*sum(log(1+R/K)) - ((K+1)/(2*K))*sum(R./(K+R));
%    h = N/(2*K*(K+1)^2) + (2*K-1)/(4*K.^2)*sum(R./(K+R)) + (K+1)/(2*K)*sum(R./(K+R).^2);
%end

end
