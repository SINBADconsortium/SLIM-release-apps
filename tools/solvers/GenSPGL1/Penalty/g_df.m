function [f,g,h] = g_df(R,K,a,flag)

N = length(R);
R = abs(R./a).^2;



if flag==1 % true
    ftrue = @(k)(-N*log(gamma((k+1)/2)/(gamma(k/2)*sqrt(pi*k))) + ((k+1)/2)*sum(log(1+R/k)));
    h = 1e-3;
    f = ftrue(K);  
    g = (ftrue(K+h) - ftrue(K-h))/(2*h);
    h = (ftrue(K+h) - 2*ftrue(K) + ftrue(K-h))/(h.^2);
    %h = N/(2*K*(K+1)^2) + (2*K-1)/(4*K.^2)*sum(R./(K+R)) + (K+1)/(2*K)*sum(R./(K+R).^2);
else % approximation
    f = (K*N/2)*log(K/(K+1)) + ((K+1)/2)*sum(log(1+R/K)) + (N/2)*(1+log(2*pi));
    g = (N/2)*(log(K/(K+1)) + 1/(K+1)) + (1/2)*sum(log(1+R/K)) - ((K+1)/(2*K))*sum(R./(K+R));
    h = N/(2*K*(K+1)^2) + (2*K-1)/(4*K.^2)*sum(R./(K+R)) + (K+1)/(2*K)*sum(R./(K+R).^2);
end
