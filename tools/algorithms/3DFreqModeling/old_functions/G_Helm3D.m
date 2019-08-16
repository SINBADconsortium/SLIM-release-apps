function M = G_Helm3D(f,m,b,o,d,n,nb,beta,l)
% Jacobian of 3D Helmholtz operator A_Helm3D
%
% d(A(m)*u)/dm = M*diag(u)
%
% use:
%   M = G_HelmD(f,m,b,o,d,n,nb,{beta},{l})
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
%   M       - sparse matrix of size n(1)*n(2)*n(3) x n(1)*n(2)*n(3)
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
% if nargin < 8
%     beta = 100;
% end
% if nargin < 9   
%     l    = 2;
% end

%% total size and angular frequency
% N  = prod(n);
w  = 2*pi*f;

%% PML
% p1 = [1 - (beta/w)*1i*linspace(1,0,nb(1)).^l    ones(1,n(1)-2*nb(1))   1 - (beta/w)*1i*linspace(0,1,nb(1)).^l]';
% p2 = [1 - (beta/w)*1i*linspace(1,0,nb(2)).^l    ones(1,n(2)-2*nb(2))   1 - (beta/w)*1i*linspace(0,1,nb(2)).^l]';
% p3 = [1 - (beta/w)*1i*linspace(1,0,nb(3)).^l    ones(1,n(3)-2*nb(3))   1 - (beta/w)*1i*linspace(0,1,nb(3)).^l]';

% p1 = repmat(permute(p1,[1 2 3]),[1 n(2) n(3)]); p1 = p1(:);
% p2 = repmat(permute(p2,[3 1 2]),[n(1) 1 n(3)]); p2 = p2(:);
% p3 = repmat(permute(p3,[2 3 1]),[n(1) n(2) 1]); p3 = p3(:);

%% Mass matrix
% combine matrices
wm1 = 0.4964958;
wm2 = 0.4510125;
wm3 = 0.052487;
wm4 = 0.45523e-5;

M = wm1*stencil([0 0 0],1,n);
M = M + (wm2/6)*stencil([1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1],[1 1 1 1 1 1],n);
M = M + (wm3/12)*stencil([1 1 0;0 1 1;1 0 1;-1 1 0;0 -1 1;-1 0 1;1 -1 0;0 1 -1;1 0 -1;-1 -1 0;0 -1 -1;-1 0 -1],[1 1 1 1 1 1 1 1 1 1 1 1],n);
M = M + (wm4/8)*stencil([1 1 1;-1 -1 -1;-1 1 1;1 -1 1;1 1 -1;-1 -1 1;1 -1 -1;-1 1 -1],[1 1 1 1 1 1 1 1],n);
a = w.^2*b;
N = length(a);
M = M*spdiags(a,0,N,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = stencil(k,a,n)

N = prod(n);
% I = @(i)sum(i(:,1)*1 + i(:,2)*n(1) + i(:,3)*n(1)*n(2),2)';

Rf = @(n)spdiags([0; ones(n-2,1); 0],0,n,n);
B  = opKron(Rf(n(3)), opKron(Rf(n(2)), Rf(n(1))));
e  = ones(N,1);

idx = sum(k(:,1)*1 + k(:,2)*n(1) + k(:,3)*n(1)*n(2),2)';

A = spalloc(N,N,length(a)*N);

if ~isempty(idx(idx==0))
    A = A + spdiags(e*a(idx==0),0,N,N);
end
if ~isempty(idx(idx~=0))
    A = A + spdiags(B*e*a(idx~=0),idx(idx~=0),N,N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function A = diags(a)
%     N = length(a);
%     A = spdiags(a,0,N,N);




