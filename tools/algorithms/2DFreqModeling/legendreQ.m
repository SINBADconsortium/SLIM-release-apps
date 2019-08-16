function [y] = legendreQ(nu,z)
% Legendre function Q_{\nu}(z) for |z| > 1 and 'small' |\nu| < O(1e2).
% 
%
% use:
%   y = legendreQ(u,nu)
%
% input:
%   u  - argument
%   nu - order
%
% output
%   y  - output, same size as u.

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

eps = 1e-2;
z(z==1) = 1-1e-10;

% use for |eta| > eps
eta = max(acosh(z),eps);
y1  = sqrt(pi)*(mygamma(1+nu)./mygamma(nu+1.5)).*exp(-(1+nu).*eta).*myhg(.5,nu+1,nu+1.5,exp(-2*eta));

% use for |eta| < eps
eta = min(acosh(z),eps);

y2 = 0*z;
a = 0.5;
b = nu + 1.5;
n = 0;
r = 1;
while(norm(r)>1e-3&&n<50)
    an = mygamma(a+n)/mygamma(a);
    bn = mygamma(b+n)/mygamma(b);
    r  = (an*bn/factorial(n).^2)*(2*mypsi(1+n) - mypsi(a+n) - mypsi(b+n)-log(1-exp(-2*eta))).*(1-exp(-2*eta)).^n;
    y2 = y2 + r;
    n  = n + 1;
end
y2 = y2.*exp(-(1+nu).*eta);

W = zeros(size(z));
W(acosh(z)>=eps) = 1;
y = W.*y1 + (1-W).*y2;
end

function y = myhg(a,b,c,z)
% hypergeometric function

y = 0*z;
r = 1;
n = 0;
while(norm(r)>1e-3&&n<50)
    an = mygamma(a+n)/mygamma(a);
    bn = mygamma(b+n)/mygamma(b);
    cn = mygamma(c+n)/mygamma(c);
    r  = (an*bn/(cn*factorial(n)))*z.^n;
    y  = y + r;
    n  = n + 1;
end
end


function [f] = mypsi(z)
%Psi     Psi (or Digamma) function valid in the entire complex plane.
%
%                 d
%        Psi(z) = --log(Gamma(z))
%                 dz
%
%usage: [f] = psi(z)
%
%tested under versions 6.0 and 5.3.1
%
%        Z may be complex and of any size.
%
%        This program uses the analytical derivative of the
%        Log of an excellent Lanczos series approximation
%        for the Gamma function.
%        
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%
%
%see also:   GAMMA GAMMALN GAMMAINC PSIN
%see also:   mhelp psi
%see also:   mhelp GAMMA

%Paul Godfrey
%pgodfrey@conexant.com
%July 13, 2001
%see gamma for calculation details...

siz = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

%reflection point
p=find(real(z)<0.5);
if ~isempty(p)
   z(p)=1-z(p);
end

%Lanczos approximation for the complex plane
 
g=607/128; % best results when 4<=g<=5
 
c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];


n=0;
d=0;
for k=size(c,1):-1:2
    dz=1./(z+k-2);
    dd=c(k).*dz;
    d=d+dd;
    n=n-dd.*dz;
end
d=d+c(1);
gg=z+g-0.5;
%log is accurate to about 13 digits...

f = log(gg) + (n./d - g./gg) ;

if ~isempty(p)
   f(p) = f(p)-pi*cot(pi*zz(p));
end

p=find(round(zz)==zz & real(zz)<=0 & imag(zz)==0);
if ~isempty(p)
   f(p) = Inf;
end

f=reshape(f,siz);

return

%A demo of this routine is:
clc
clear all
close all
x=-4:1/16:4.5;
y=-4:1/16:4;
[X,Y]=meshgrid(x,y);
z=X+i*Y;
f=psi(z);
p=find(abs(f)>10);
f(p)=10;

mesh(x,y,abs(f),phase(f));
view([45 10]);
rotate3d;

figure(2);
ezplot psi;
grid on;

One=psi(2)-psi(1)
EulerGamma=-psi(1)

end
function [f] = mygamma(z)
% GAMMA  Gamma function valid in the entire complex plane.
%        Accuracy is 15 significant digits along the real axis
%        and 13 significant digits elsewhere.
%        This routine uses a superb Lanczos series
%        approximation for the complex Gamma function.
%
%        z may be complex and of any size.
%        Also  n! = prod(1:n) = gamma(n+1)
%
%usage: [f] = gamma(z)
%       
%tested on versions 6.0 and 5.3.1 under Sun Solaris 5.5.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931-944
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%            W. J. Cody "An Overview of Software Development for Special
%            Functions", 1975
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%
%Paul Godfrey
%pgodfrey@intersil.com
%http://winnie.fit.edu/~gabdo/gamma.txt
%Sept 11, 2001

siz = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

p=find(real(z)<0);
if ~isempty(p)
   z(p)=-z(p);
end

% 15 sig. digits for 0<=real(z)<=171
% coeffs should sum to about g*g/2+23/24

g=607/128; % best results when 4<=g<=5

c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];

%Num Recipes used g=5 with 7 terms
%for a less effective approximation

z=z-1;
zh =z+0.5;
zgh=zh+g;
%trick for avoiding FP overflow above z=141
zp=zgh.^(zh*0.5);

ss=0.0;
for pp=size(c,1)-1:-1:1
    ss=ss+c(pp+1)./(z+pp);
end

%sqrt(2Pi)
sq2pi=  2.5066282746310005024157652848110;
f=(sq2pi*(c(1)+ss)).*((zp.*exp(-zgh)).*zp);

f(z==0 | z==1) = 1.0;

%adjust for negative real parts
if ~isempty(p)
   f(p)=-pi./(zz(p).*f(p).*sin(pi*zz(p)));
end

%adjust for negative poles
p=find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
   f(p)=Inf;
end

f=reshape(f,siz);

end