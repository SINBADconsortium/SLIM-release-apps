function [H,dH_forw,dH_adj,DdH_adj] = Poiss2D(epsilon,n)
% Finite volume discretization of the Poisson equation
%
%   div(epsilon grad(u)) = q
% 
%   
% Curt Da Silva, 2016
%
% Usage:
%   [H,dH_forw,dH_adj] = Poiss2D(epsilon,n);
%
% Input:
%   epsilon - (nz+2) x (nx+2) 
% 
% Output:
%   H       - nz*nx x nz*nx system matrix
%   dH_forw - function of (deps,u), derivative mapping deps -> DH(eps)[deps]u
%   dH_adj  - function of (z,u), adjoint of dH_forw
    

    assert(numel(epsilon)==prod(n(1:2)),'epsilon must be nz x nx');
    nz = n(1); nx = n(2); 
    
    int_z = ones(nz,1); int_x = ones(nx,1);
    %int_z([1 end]) = 0; int_x([1 end]) = 0;
    [IZ,IX] = ndgrid(int_z,int_x);
    
    Rext = opKron(opExtension(nx,[1 1],1),opExtension(nz,[1 1],1));
    Rext_adj = opKron(opRestriction(nx+2,2:(nx+1)),opRestriction(nz+2,2:(nz+1)));
    e = Rext*vec(epsilon);
    Rz = @(offset) opRestriction(nz+2,(2+offset):(nz+1+offset));
    Rx = @(offset) opRestriction(nx+2,(2+offset):(nx+1+offset));
    R = @(oz,ox) opKron(Rx(ox),Rz(oz));        
    
    b0 = -(R(0,0)+R(-1,0)+R(0,-1)+R(-1,-1));
    b1 = opDiag_swp(vec(IZ))*(1/2*(R(0,0)+R(0,-1)));
    b2 = opDiag_swp(vec(IX))*(1/2*(R(0,0)+R(-1,0)));
    b3 = opDiag_swp(vec(IZ))*(1/2*(R(-1,-1)+R(-1,0)));
    b4 = opDiag_swp(vec(IX))*(1/2*(R(0,-1)+R(-1,-1)));
    e1 = ones(prod(n),1);    
    
    A = @(oz) spdiags(e1,oz,prod(n),prod(n));
    sdiag = @(x) spdiags(x,0,length(x),length(x));
    H = A(0)*sdiag(b0*e) + A(1)*sdiag(b1*e) + A(nz)*sdiag(b2*e) + A(-1)*sdiag(b3*e) + A(-nz)*sdiag(b4*e);
    dH_forw = @(de,u) (A(0)*sdiag(b0*(Rext*de)) + A(1)*sdiag(b1*(Rext*de)) + A(nz)*sdiag(b2*(Rext*de)) + A(-1)*sdiag(b3*(Rext*de)) + A(-nz)*sdiag(b4*(Rext*de)))*u;
    dH_adj = @(z,u) Rext_adj*(b0'*(conj(u) .* (A(0)'*z)) + b1'*(conj(u) .* (A(1)'*z)) + b2'*(conj(u) .* (A(nz)'*z)) + b3'*(conj(u) .* (A(-1)'*z)) + b4'*(conj(u) .* (A(-nz)'*z)));
    DdH_adj = @(u,du,dm,z) Rext_adj*(b0'*(conj(du) .* (A(0)'*z)) + b1'*(conj(du) .* (A(1)'*z)) + b2'*(conj(du) .* (A(nz)'*z)) + b3'*(conj(du) .* (A(-1)'*z)) + b4'*(conj(du) .* (A(-nz)'*z))); 
    
end

function [H,dH_forw,dH_adj] = Poiss2D_fd(epsilon,n)
    nz = n(1); nx = n(2); N = nz*nx;
    % indicator vectors for the interior nodes
    int_z = ones(nz,1); int_x = ones(nx,1);
    int_z([1 end]) = 0; int_x([1 end]) = 0;
    [IZ,IX] = ndgrid(int_z,int_x);
    
    % Extend epsilon by one grid point in each direction
    Rext = opKron(opExtension(nx,[1 1],1),opExtension(nz,[1 1],1));
    eps_ext = Rext*vec(epsilon);
    
    Rz = @(offset) opRestriction_swp(nz+2,(2+offset):(nz+1+offset));
    Rx = @(offset) opRestriction_swp(nx+2,(2+offset):(nx+1+offset));
    R = @(oz,ox) opKron(Rx(ox),Rz(oz)); 
    aPN = vec(IZ) .* (R(1,0)*eps_ext);
    aNN = -(R(1,0)*eps_ext + R(-1,0)*eps_ext ) - (R(0,1)*eps_ext + R(0,-1)*eps_ext);
    aMN = vec(IZ) .* (R(-1,0)*eps_ext);
    aNP = vec(IX) .* (R(0,1)*eps_ext);
    aNM = vec(IX) .* (R(0,-1)*eps_ext);
    
    H = spdiags([aPN,aMN,aNN,aNP,aNM],[1,-1,0,nz,-nz],N,N);
    
    dH_forw = @(de,u) dH_func(de,u,R,Rext,IZ,IX);
    
    dH_adj = @(z,u) dH_adj_func(z,u,R,Rext,IZ,IX);
end

function dH = dH_func(de,u,R,Rext,IZ,IX)
    nz = size(IZ,1); nx = size(IZ,2);
    dex = vec(Rext*de);
    aPN = vec(IZ) .* (R(1,0)*dex);
    aNN = -(R(1,0)*dex + R(-1,0)*dex ) - (R(0,1)*dex + R(0,-1)*dex);
    aMN = vec(IZ) .* (R(-1,0)*dex );
    aNP = vec(IX) .* (R(0,1)*dex);
    aNM = vec(IX) .* (R(0,-1)*dex);
    
    dH = (spdiags([aPN,aMN,aNN,aNP,aNM],[1,-1,0,nz,-nz],nz*nx,nz*nx))*u;
end

function y = dH_adj_func(z,u,R,Rext,IZ,IX)
    p = size(u,2);
    nz = size(IZ,1); 
    % aPN        
    y = Rext'*(R(1,0)'*(repmat(vec(IZ),1,p).* (conj(u).*circshift(z,1))));
    % aMN
    y = y + Rext'*(R(-1,0)'*(repmat(vec(IZ),1,p).* (conj(u).*circshift(z,-1))));
    % aNN
    y = y - Rext'*(R(1,0)+R(-1,0)+R(0,1)+R(0,-1))'*(conj(u) .* z);
    % aNP
    y = y + Rext'*(R(0,1)'*(repmat(vec(IX),1,p).* (conj(u).*circshift(z,nz))));
    % aNM
    y = y + Rext'*(R(0,-1)'*(repmat(vec(IX),1,p).* (conj(u).*circshift(z,-nz))));
end
