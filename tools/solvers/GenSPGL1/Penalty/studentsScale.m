function [y g h] = studentsScale(sig2, para)
% function [y g h] = students(w, para)

% Input:  sig2 = sigma^2 parameter
%         para.n  number of measurements 
%         para.r  residual vector (length n)


n = para.n;
r = para.r;
k = para.k;


y = (n/2)*log(sig2) + 0.5*(k+1)*sum(log(1 + r.^2./(k*sig2)));

%g = -n*(k)/(2*sig2) + 0.5*k*(k+1)*sum(1./(r.^2+k*sig2)); 
g = n/(2*sig2) - ((k+1)/(2*k*sig2))*sum(r.^2./(sig2 + r.^2/k)); 

%h = n*(k)/(2*sig2^2) - 0.5*k^2*(k+1)*sum(1./(r.^2+k*sig2).^2);
h = -n/(2*sig2^2) + ((k+1)/(2*k*sig2^2))*sum(r.^2./(sig2 + r.^2/k)) +  ((k+1)/(2*k*sig2))*sum(r.^2./(sig2 + r.^2/k).^2);


end
