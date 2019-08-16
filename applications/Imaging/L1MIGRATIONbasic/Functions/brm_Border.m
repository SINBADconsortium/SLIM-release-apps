function  [b,nsr] = brm_Border(xd,zd,nx,nz,nx_border,nz_border,ns,p,israndn) 
% function  [b,nsr] = brm_Border(xd,zd,nx,nz,nx_border,nz_border,ns,p,israndn) 
%  
% this opertor sets simutanous source function b,
% p is your sampeling rate
% israndn: if 0, put sources periodiclly; if 1, randomly put sources.
% 
if ~exist('israndn'),israndn = 0;end




hx = xd/(nx-1);
hz = zd/(nz-1);
% Q  = zeros(nz,nx,ns);
% Q(1,1:nx,1:ns) = eye(ns)/(hx*hz);
Q  = spalloc(nz,nx*ns,nx);

if israndn
    idx = randperm(nx);
    idx = sort(idx(1:ns));
    idx = sub2ind([ns,nx],1:ns,idx);
    idx  = idx(:);
else
    n  = (nx-1)./ns;
    idx= round(1:nx+n:nx*ns);
end

Q(1,idx) = 1/(hx*hz);
Q  = Q(:);


m = nx;
n = ns;
%p = 0.80;
k = round(p * n);

% define operator
%mask                 = (randperm(n) > k);
%Rs                   = opColumnRestrict(m,n,find(mask),'zeros');
%Ms                   = opKron(opGaussian(n,n),opDirac(m));
%Iz                   = opDirac(nz);
%Ix                   = opDirac(nx);
%As                   = opFoG(Rs,Ms);
%A                    = opKron(As,Iz);


% take measurements
%RMa = A(Q(:),1);
%RMa = reshape(RMa,nz,nx,ns);

%figure(1);
%imagesc(reshape(Q(:,:,50),nz,nx));colormap('gray')
%figure(2);
%imagesc(reshape(RMa(:,:,50),nz,nx));colormap('gray')

% define operator
mask                 = (randperm(n) > k);
Rs                   = opColumnRestrict(m,n,find(mask),'discard');
Ms                   = opKron(opGaussian(n,n),opDirac(m));
Iz                   = opDirac(nz);
Ix                   = opDirac(nx);
As                   = opFoG(Rs,Ms);
A                    = opKron(As,Iz);

% take measurements
RMa = A*Q;
sa  = size(Rs);
nsr = sa(1)/nx;
RMa = reshape(RMa,nz,nx,nsr);

% add border to the arguement


RMa_border           = zeros(nz+2*nz_border,nx+2*nx_border,nsr);
nsrloop              = 1:nsr;
RMa_border(nz_border+1:nz_border+nz,nx_border+1:nx_border+nx,nsrloop) = RMa; 
b = spalloc((nx+2*nx_border)*(nz+2*nz_border),nsr,nx*nsr);
% b(1:nx,1:nsr) = RMa(1,:,:);
b = reshape(RMa_border,(nx+2*nx_border)*(nz+2*nz_border),nsr);

whos b
