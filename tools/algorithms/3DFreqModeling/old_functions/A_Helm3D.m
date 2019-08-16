function [A] = A_Helm3D(f,m,b,o,d,n,nb,beta,l)
% 3D Helmholtz operator, discretized using 27-point stencil with PML
% [Operto at al. 2007. Geophysics 72(5), SM195]
%
% The PML functions are of the form: xi = 1 + i*(beta/w)*[0:1/(nb-1):1].^l
%
% Tristan van Leeuwen, 2012
% tleeuwen@eos.ubc.ca
%
% use:
%   A = A_HelmD(f,m,b,o,d,n,nb,{beta},{l})
%
% input:
%   f       - frequency [1/s]
%   m       - gridded slowness-squared [s^2/m^2]
%   b       - gridded buoyancy 1/rho [cm^3/gr]
%   {o,d,n} - grid definition: z = o(1) + [0:n(1)-1]*d(1), etc.
%   nb      - number of points to use as absorbing boundary in each direction
%   {beta,l}- PML parameters, default beta = 100, l = 2.
%
% output:
%   A       - sparse matrix of size n(1)*n(2)*n(3) x n(1)*n(2)*n(3)
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

% input checking
if nargin < 8
    beta = 100;
end
if nargin < 9    
    l    = 2;
end

%% total size and angular frequency
N  = prod(n);
w  = 2*pi*f;

%% PML
p1 = [1 - (beta/w)*1i*linspace(1,0,nb(1)).^l    ones(1,n(1)-2*nb(1))   1 - (beta/w)*1i*linspace(0,1,nb(1)).^l]';
p2 = [1 - (beta/w)*1i*linspace(1,0,nb(2)).^l    ones(1,n(2)-2*nb(2))   1 - (beta/w)*1i*linspace(0,1,nb(2)).^l]';
p3 = [1 - (beta/w)*1i*linspace(1,0,nb(3)).^l    ones(1,n(3)-2*nb(3))   1 - (beta/w)*1i*linspace(0,1,nb(3)).^l]';

p1 = repmat(permute(p1,[1 2 3]),[1 n(2) n(3)]); p1 = p1(:);
p2 = repmat(permute(p2,[3 1 2]),[n(1) 1 n(3)]); p2 = p2(:);
p3 = repmat(permute(p3,[2 3 1]),[n(1) n(2) 1]); p3 = p3(:);

%% averaging matrices
%Av = @(k)kron(Af(k(3),n(3)), kron(Af(k(2),n(2)),Af(k(1),n(1))));
Av = @(k)opKron(Af(k(3),n(3)), Af(k(2),n(2)),Af(k(1),n(1)));

%% Mass matrix
% combine matrices
wm1 = 0.4964958;
wm2 = 0.4510125;
wm3 = 0.052487;
wm4 = 0.45523e-5;

M = wm1*stencil([0 0 0],1,n)...
 + (wm2/6)*stencil([1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1],[1 1 1 1 1 1],n)...
 + (wm3/12)*stencil([1 1 0;0 1 1;1 0 1;-1 1 0;0 -1 1;-1 0 1;1 -1 0;0 1 -1;1 0 -1;-1 -1 0;0 -1 -1;-1 0 -1],[1 1 1 1 1 1 1 1 1 1 1 1],n)...
 + (wm4/8)*stencil([1 1 1;-1 -1 -1;-1 1 1;1 -1 1;1 1 -1;-1 -1 1;1 -1 -1;-1 1 -1],[1 1 1 1 1 1 1 1],n);

%% Stifness matrix
bx = @(k)(Av(k)*b)./(Av(k)*p1)./(d(1).^2*p1);
by = @(k)(Av(k)*b)./(Av(k)*p2)./(d(2).^2*p2);
bz = @(k)(Av(k)*b)./(Av(k)*p3)./(d(3).^2*p3);

% Sc,  eq. (D-2)
Ix = [1 0 0;-1  0  0]; 
Iy = [0 1 0; 0 -1  0];
Iz = [0 0 1; 0  0 -1];
Sc = getS(@(k)bx(k),Ix,n) +  getS(@(k)by(k),Iy,n) + getS(@(k)bz(k),Iz,n);
% Sx
Ix = [ 1  0  0;-1  0  0]; 
Iy = [ 0  1  1; 0 -1 -1; 0  1 -1; 0 -1  1];
Iz = [ 0  1  1; 0 -1 -1; 0 -1  1; 0  1 -1];
Sx =  getS(@(k)bx(k),Ix,n) + getS(@(k).25*by(k),Iy,n) + getS(@(k).25*bz(k),Iz,n);
% Sy
Ix = [ 1  0  1;-1  0 -1; 1 0 -1;-1 0  1]; 
Iy = [ 0  1  0; 0 -1  0];
Iz = [ 1  0  1;-1  0 -1;-1 0  1; 1 0 -1];
Sy =  getS(@(k).25*bx(k),Ix,n) + getS(@(k)by(k),Iy,n) + getS(@(k).25*bz(k),Iz,n); 
% Sz
Ix = [ 1  1  0;-1 -1  0; 1 -1  0;-1  1  0]; 
Iy = [ 1  1  0;-1 -1  0;-1  1  0; 1 -1  0];
Iz = [ 0  0  1; 0  0 -1];
Sz =  getS(@(k).25*bx(k),Ix,n) + getS(@(k).25*by(k),Iy,n) + getS(@(k)bz(k),Iz,n);
% S1
Ix = [ 1  1 -1;-1 -1  1; 1 -1  1;-1  1 -1]; 
Iy = [ 1  1  1;-1 -1 -1;-1  1 -1; 1 -1  1];
Iz = [ 1  1  1;-1 -1 -1;-1 -1  1; 1  1 -1];
S1 =  getS(@(k).25*bx(k),Ix,n) + getS(@(k).25*by(k),Iy,n) + getS(@(k).25*bz(k),Iz,n);
% S2
Ix = [ 1  1  1;-1 -1 -1; 1 -1 -1;-1  1  1]; 
Iy = [ 1  1 -1;-1 -1  1;-1  1  1; 1 -1 -1];
Iz = [ 1  1  1;-1 -1 -1;-1 -1  1; 1  1 -1];
S2 =  getS(@(k).25*bx(k),Ix,n) + getS(@(k).25*by(k),Iy,n) + getS(@(k).25*bz(k),Iz,n); 
% S3
Ix = [ 1  1  1;-1 -1 -1; 1 -1 -1;-1  1  1]; 
Iy = [ 1  1  1;-1 -1 -1;-1  1 -1; 1 -1  1];
Iz = [ 1 -1  1;-1  1 -1;-1  1  1; 1 -1 -1];
S3 =  getS(@(k).25*bx(k),Ix,n) + getS(@(k).25*by(k),Iy,n) + getS(@(k).25*bz(k),Iz,n);
% S4
Ix = [ 1  1 -1;-1 -1  1; 1 -1  1;-1  1 -1]; 
Iy = [ 1  1 -1;-1 -1  1;-1  1  1; 1 -1 -1];
Iz = [ 1 -1  1;-1  1 -1;-1  1  1; 1 -1 -1];
S4 =  getS(@(k).25*bx(k),Ix,n) + getS(@(k).25*by(k),Iy,n) + getS(@(k).25*bz(k),Iz,n);
%%
% Combine matrices
w1 = 1.8395265e-5;
w2 = 0.890077;
w3 = 0.1099046;

%% final matrix
A = M*diags(w^2*m.*b) + w1*Sc + (w2/3)*(Sx + Sy + Sz) + (w3/4)*(S1 + S2 + S3 + S4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Af(k,n)
if k==1
	A = spdiags(ones(n,1)*[.5 .5],[0 1],n,n);
elseif k==-1
	A = spdiags(ones(n,1)*[.5 .5],[-1 0],n,n);
else
	A = spdiags(ones(n,1),0,n,n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Df(k,n)
if k==1
	A = spdiags(ones(n,1)*[-1 1],[0 1],n,n);
elseif k==-1
	A = spdiags(ones(n,1)*[-1 1],[-1 0],n,n);
else
	A = spdiags(ones(n,1),0,n,n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = stencil(k,a,n)

N = prod(n);
I = @(i)sum(i(:,1)*1 + i(:,2)*n(1) + i(:,3)*n(1)*n(2),2)';

Rf = @(n)spdiags([0; ones(n-2,1); 0],0,n,n);
B  = kron(Rf(n(3)), kron(Rf(n(2)), Rf(n(1))));
e  = ones(N,1);

idx = I(k);

A = spalloc(N,N,length(a)*N);

if ~isempty(idx(idx==0))
    A = A + spdiags(e*a(idx==0),0,N,N);
end
if ~isempty(idx(idx~=0))
    A = A + spdiags(B*e*a(idx~=0),idx(idx~=0),N,N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = diags(a)
    N = length(a);
    A = spdiags(a,0,N,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = getS(b,I,n)

m = size(I,1);
w = sign(mod(1:4,2)-.5);
k = 1;

Ik = I/2 + ones(m,1)*I(k,:)/2;
S  =  w(k)*diags(b(I(k,:)))*stencil(Ik,w,n);
    
for k = 2:m
    Ik = I/2 + ones(m,1)*I(k,:)/2;
    S = S + w(k)*diags(b(I(k,:)))*stencil(Ik,w,n);
end
    



