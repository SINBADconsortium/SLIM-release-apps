function [R,idx] = R_Helm3D(f,m,b,o,d,n,nb,beta,l)
% 3D Helmholtz operator, discretized using 27-point stencil with PML
% [Operto at al. 2007. Geophysics 72(5), SM195]
%
% The PML functions are of the form: xi = 1 + i*(beta/w)*[0:1/(nb-1):1].^l
% use:
%   [R,idx] = R_HelmD(f,k,b,o,d,n,nb,{beta},{l})
%
% input:
%   f       - frequency [1/s]
%   m       - gridded squared-slowness [s^2/m^2]
%   b       - gridded buoyancy 1/rho [cm^3/gr]
%   {o,d,n} - grid definition: z = o(1) + [0:n(1)-1]*d(1), etc.
%   nb      - number of points to use as absorbing boundary in each direction
%   {beta,l}- PML parameters, default beta = 100, l = 2.
%
% output:
%   R       - matrix in band-storate format
%   idx     - indices
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
%p3 = [ones(1,n(3))]';

p1 = repmat(permute(p1,[1 2 3]),[1 n(2) n(3)]); p1 = p1(:);
p2 = repmat(permute(p2,[3 1 2]),[n(1) 1 n(3)]); p2 = p2(:);
p3 = repmat(permute(p3,[2 3 1]),[n(1) n(2) 1]); p3 = p3(:);

%% averaging matrices
Av = @(k)opKron(Af(k(3),n(3)), Af(k(2),n(2)),Af(k(1),n(1)));

%% allocote memory
R = zeros(N,27);

%% Mass matrix
wm1 = 0.4964958;
wm2 = 0.4510125;
wm3 = 0.052487;
wm4 = 0.45523e-5;

R = stencil_add(R,[0 0 0],wm1,w^2*m.*b);
I = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
R = stencil_add(R,I,(wm2/6),shift(w^2*m.*b,I,n));
I = [1 1 0;0 1 1;1 0 1;-1 1 0;0 -1 1;-1 0 1;1 -1 0;0 1 -1;1 0 -1;-1 -1 0;0 -1 -1;-1 0 -1];
R = stencil_add(R,I,(wm3/12),shift(w^2*m.*b,I,n));
I = [1 1 1;-1 -1 -1;-1 1 1;1 -1 1;1 1 -1;-1 -1 1;1 -1 -1;-1 1 -1];
R = stencil_add(R,I,(wm4/8),shift(w^2*m.*b,I,n));

%% Stifness matrix
w1 = 1.8395265e-5;
w2 = 0.890077;
w3 = 0.1099046;

bx = @(k)(Av(k)*b)./(Av(k)*p1)./(d(1).^2*p1);
by = @(k)(Av(k)*b)./(Av(k)*p2)./(d(2).^2*p2);
bz = @(k)(Av(k)*b)./(Av(k)*p3)./(d(3).^2*p3);

% Sc,  eq. (D-2)
Ix = [1 0 0;-1  0  0]; 
Iy = [0 1 0; 0 -1  0];
Iz = [0 0 1; 0  0 -1];
R = getS(R,@(k)w1*bx(k),Ix,n); 
R = getS(R,@(k)w1*by(k),Iy,n);
R = getS(R,@(k)w1*bz(k),Iz,n);
% Sx
Ix = [ 1  0  0;-1  0  0]; 
Iy = [ 0  1  1; 0 -1 -1; 0  1 -1; 0 -1  1];
Iz = [ 0  1  1; 0 -1 -1; 0 -1  1; 0  1 -1];
R =  getS(R,@(k)(w2/3)*bx(k),Ix,n);
R =  getS(R,@(k).25*(w2/3)*by(k),Iy,n);
R =  getS(R,@(k).25*(w2/3)*bz(k),Iz,n);
% Sy
Ix = [ 1  0  1;-1  0 -1; 1 0 -1;-1 0  1]; 
Iy = [ 0  1  0; 0 -1  0];
Iz = [ 1  0  1;-1  0 -1;-1 0  1; 1 0 -1];
R = getS(R,@(k).25*(w2/3)*bx(k),Ix,n);
R = getS(R,@(k)(w2/3)*by(k),Iy,n);
R = getS(R,@(k).25*(w2/3)*bz(k),Iz,n); 
% Sz
Ix = [ 1  1  0;-1 -1  0; 1 -1  0;-1  1  0]; 
Iy = [ 1  1  0;-1 -1  0;-1  1  0; 1 -1  0];
Iz = [ 0  0  1; 0  0 -1];
R = getS(R,@(k).25*(w2/3)*bx(k),Ix,n);
R = getS(R,@(k).25*(w2/3)*by(k),Iy,n);
R = getS(R,@(k)(w2/3)*bz(k),Iz,n);
% S1
Ix = [ 1  1 -1;-1 -1  1; 1 -1  1;-1  1 -1]; 
Iy = [ 1  1  1;-1 -1 -1;-1  1 -1; 1 -1  1];
Iz = [ 1  1  1;-1 -1 -1;-1 -1  1; 1  1 -1];
R = getS(R,@(k).25*(w3/4)*bx(k),Ix,n);
R = getS(R,@(k).25*(w3/4)*by(k),Iy,n);
R = getS(R,@(k).25*(w3/4)*bz(k),Iz,n);
% S2
Ix = [ 1  1  1;-1 -1 -1; 1 -1 -1;-1  1  1]; 
Iy = [ 1  1 -1;-1 -1  1;-1  1  1; 1 -1 -1];
Iz = [ 1  1  1;-1 -1 -1;-1 -1  1; 1  1 -1];
R = getS(R,@(k).25*(w3/4)*bx(k),Ix,n);
R = getS(R,@(k).25*(w3/4)*by(k),Iy,n);
R = getS(R,@(k).25*(w3/4)*bz(k),Iz,n);
% S3
Ix = [ 1  1  1;-1 -1 -1; 1 -1 -1;-1  1  1]; 
Iy = [ 1  1  1;-1 -1 -1;-1  1 -1; 1 -1  1];
Iz = [ 1 -1  1;-1  1 -1;-1  1  1; 1 -1 -1];
R = getS(R,@(k).25*(w3/4)*bx(k),Ix,n);
R = getS(R,@(k).25*(w3/4)*by(k),Iy,n);
R = getS(R,@(k).25*(w3/4)*bz(k),Iz,n);
% S4
Ix = [ 1  1 -1;-1 -1  1; 1 -1  1;-1  1 -1]; 
Iy = [ 1  1 -1;-1 -1  1;-1  1  1; 1 -1 -1];
Iz = [ 1 -1  1;-1  1 -1;-1  1  1; 1 -1 -1];
R = getS(R,@(k).25*(w3/4)*bx(k),Ix,n);
R = getS(R,@(k).25*(w3/4)*by(k),Iy,n);
R = getS(R,@(k).25*(w3/4)*bz(k),Iz,n);
%%


Rf = @(n)spdiags([0; ones(n-2,1); 0],0,n,n);
B  = kron(Rf(n(3)), kron(Rf(n(2)), Rf(n(1))));
e  = B*ones(N,1);
I = find(e==0);

idx = [];
for k=[-1:1],for l=[-1:1],for m=[-1:1] idx = [idx m + n(1)*l + n(1)*n(2)*k];end,end,end

for k=[1:13 15:27]
    Ik = mod(I-1-idx(k),N)+1;
    R(Ik,k) = 0;    
end
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
function S = getS(S,b,I,n)

N = prod(n);
m = size(I,1);
w = sign(mod(1:m,2)-.5);
    
for k = 1:m
    Ik = I/2 + ones(m,1)*I(k,:)/2;
    ak = w(k)*b(I(k,:));
    S = stencil_add(S,Ik,ak,w);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = stencil_add(S,I,a,w)
i = idx(I);
S(:,i) = S(:,i) + a*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = idx(I)
i = sum(I(:,1)*1 + I(:,2)*3 + I(:,3)*9,2)'+14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = shift(e,I,n)

m = size(I,1);
S = zeros(prod(n),m);

for k = 1:m
    i = sum(I(k,:).*[1 n(1) n(1)*n(2)]);
    if i > 0
        S(:,k) = [e(i+1:end);zeros(i,1)]; 
    else
        i = abs(i);
        S(:,k) = [zeros(i,1);e(1:end-i)];
    end
end
