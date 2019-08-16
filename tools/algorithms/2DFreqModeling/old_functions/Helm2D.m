function [H,A,Bk,C] = Helm2D(k,o,d,n,nb)
% 2D Helmholtz operator, discretized using `optimal' 9-point stencil [Jo et al '96, Geophysics 62(2), 529-537]
% and damping boundary layer in all directions
%
% WORKS ONLY FOR SQUARE GRIDS
%
% use:
%   [H,A,B] = Helm2D(k,o,d,n,nb)
%
% input:
%   k       - gridded wavenumber 2*pi*f/c [1/m]
%   {o,d,n} - grid definition: z = o(1) + [0:n(1)-1]*d(1), etc.
%   nb      - number of points to use as absorbing boundary in each direction
%
% output:
%   H       - Helmholtz matrix of size n(1)*n(2) x n(1)*n(2)
%   A,B     - H = A * diag( (B*k).^2 ) + const, matrices used in derivative calculation
%
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth & Ocean Sciences
%         The University of British Columbia
% 
% Date: February, 2012
%
% Updated by : Curt Da Silva
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth & Ocean Sciences
%         The University of British Columbia
%
% Date : January, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

N = prod(n);
beta = .15;

% boundary layer
w1 = ones(n(1),1);
w2 = ones(n(2),1);
if nb(1)
    r1 = [linspace(1,0,nb(1))';zeros(n(1)-2*nb(1),1);linspace(0,1,nb(1))'];
    w1 = 1./(1+1i*beta*r1.^2);
end
if nb(2)
    r2 = [linspace(1,0,nb(2))';zeros(n(2)-2*nb(2),1);linspace(0,1,nb(2))'];
    w2 = 1./(1+1i*beta*r2.^2);
end

Bk = kron(spdiags(w2,0,n(2),n(2)),spdiags(w1,0,n(1),n(1)));
k  = Bk*k;

% 5-point stencil coeff.
c1 = ones(N,1)*[1 -2 1]/d(1)^2; % d^2/dz^2 
c2 = ones(N,1)*[1 -2 1]/d(2)^2; % d^2/dx^2

% rotated 5-point stencil coeff.
h   = sqrt(d(1)^2 + d(2)^2);
rc1 = ones(N,1)*[1 -2 1]/h^2;
rc2 = ones(N,1)*[1 -2 1]/h^2;

% BC's
% We use Sommerfield radiation conditions for $L^(0)$ and set the L^{(45)}
% coeff. to zero.

%define boundary indices
%top
Ia1 = 1:n(1):N;
%bottom
Ib1 = n(1):n(1):N;
%left
Ia2 = 1:n(1);
%right
Ib2 = (n(2)-1)*n(1) + [1:n(1)];

% top
c1(Ia1,3) = 0; 
c1(Ia1,2) = c1(Ia1,2) + (1./(1+1i*k(Ia1')*d(1)))/d(1)^2;
rc1(Ia1,:) = 0;
rc2(Ia1,:) = 0;

% bottom
c1(Ib1,1) = 0; 
c1(Ib1,2) = c1(Ib1,2) + (1./(1+1i*k(Ib1')*d(1)))/d(1)^2;
rc1(Ib1,:) = 0;
rc2(Ib1,:) = 0;

% left
c2(Ia2,2) = c2(Ia2,2) + (1./(1+1i*k(Ia2')*d(2)))/d(2)^2;
rc1(Ia2,:) = 0;
rc2(Ia2,:) = 0;

% right
c2(Ib2,2) = c2(Ib2,2) + (1./(1+1i*k(Ib2')*d(2)))/d(2)^2;
rc1(Ib2,:) = 0;
rc2(Ib2,:) = 0;


% assemble coeff. into matrices
L0  = spdiags(c1,[-1 0 1],N,N) + spdiags(c2,[-n(1) 0 n(1)],N,N);
L45 = spdiags(rc1,[-n(1)+1 0 n(1)-1],N,N) + spdiags(rc2,[-n(1)-1 0 n(1)+1],N,N);
A0  = spdiags(ones(N,4),[-n(1) -1 1 n(1)],N,N);
A45 = spdiags(ones(N,4),[-n(1)-1 -n(1)+1 n(1)-1 n(1)+1],N,N);
E   = spdiags(k.^2,0,N,N);

% define params
alpha = 0.5461; delta = 0.09381; gamma = 0.6248; % optimized according to Jo et al '96, Geophysics 62(2), 529-537
%gamma = 1; alpha = 1; delta = 0; %gives normal 5-point discretization

A = (gamma*speye(N) + delta*A0 + .25*(1-gamma-4*delta)*A45);
C = alpha*L0 + (1-alpha)*L45;
H  = A*E + C;

