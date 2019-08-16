function [RR,idx] = pR_Helm3D(f,m,b,o,d,n,nb,beta,l)
% Distributed 3D Helmholtz operator, discretized using 27-point stencil with PML
% [Operto at al. 2007. Geophysics 72(5), SM195]
%
% The PML functions are of the form: xi = 1 + i*(beta/w)*[0:1/(nb-1):1].^l
%
% use:
%   [RR,idx] = pR_HelmD(f,m,b,o,d,n,nb,{beta},{l})
%
% input:
%   f       - frequency [1/s]
%   m       - gridded squared-slowness [s^2/m^2], it should be a distributed
%             3D cube, not an array. It should be distributed in the last dimesion.
%   b       - gridded buoyancy 1/rho [cm^3/gr], it should be a distributed
%             3D cube, not an array. It should be distributed in the last dimesion.
%   {o,d,n} - grid definition: z = o(1) + [0:n(1)-1]*d(1), etc.
%   nb      - number of points to use as absorbing boundary in each direction
%   {beta,l}- PML parameters, default beta = 100, l = 2.
%
% output:
%   RR      - matrix in band-storate format, distributed
%   idx     - indices
%
% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: April, 2014
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
w  = 2*pi*f;

spmd
    %  For P
    id    = labindex;
    nlab  = numlabs;
    ml    = getLocalPart(m);
    bl    = getLocalPart(b);
    if id < nlab
        labSend(ml(:,:,end-1:end),id+1);
    end
    if id > 1
        m1 = labReceive(id-1);
        labSend(ml(:,:,1:2),id-1);
    end
    if id < nlab
        me = labReceive(id+1);
    end
    
    % Set the local size
    codist    = getCodistributor(m);
    Part      = codist.Partition;
    nl        = [n(1), n(2), Part(id)];
    
    % Set the local model
    if id == 1
        mln 			 = zeros(nl(1),nl(2),nl(3)+2);
    else if id == nlab
            mln              = zeros(nl(1),nl(2),nl(3)+2);
        else
            mln 			 = zeros(nl(1),nl(2),nl(3)+4);
        end
    end
    
    
    %%%%%%%%
    if id < nlab
        labSend(bl(:,:,end-1:end),id+1);
    end
    if id > 1
        b1 = labReceive(id-1);
        labSend(bl(:,:,1:2),id-1);
    end
    if id < nlab
        be = labReceive(id+1);
    end
    
    
    % Set the local model
    if id == 1
        bln 			 = zeros(nl(1),nl(2),nl(3)+2);
    else if id == nlab
            bln              = zeros(nl(1),nl(2),nl(3)+2);
        else
            bln 			 = zeros(nl(1),nl(2),nl(3)+4);
        end
    end
    
    
    if id == 1
        mln(:,:,1:end-2)       = ml;
        mln(:,:,end-1:end)     = me;
        bln(:,:,1:end-2)       = bl;
        bln(:,:,end-1:end)     = be;
        nl(3)                  = nl(3)+2;
    else if id == nlab
            mln(:,:,1:2)       = m1;
            mln(:,:,3:end)     = ml;
            bln(:,:,1:2)       = b1;
            bln(:,:,3:end)     = bl;
            nl(3)              = nl(3)+2;
        else
            mln(:,:,3:end-2)       = ml;
            mln(:,:,1:2)           = m1;
            mln(:,:,end-1:end)     = me;
            bln(:,:,3:end-2)       = bl;
            bln(:,:,1:2)           = b1;
            bln(:,:,end-1:end)     = be;
            nl(3)                  = nl(3)+4;
        end
    end
    
    bl = [];
    ml = [];
    
    mln = mln(:);
    bln = bln(:);
  
    N                = prod(n);
    
    %% PML
    p1 = [1 - (beta/w)*1i*linspace(1,0,nb(1)).^l    ones(1,n(1)-2*nb(1))   1 - (beta/w)*1i*linspace(0,1,nb(1)).^l]';
    p2 = [1 - (beta/w)*1i*linspace(1,0,nb(2)).^l    ones(1,n(2)-2*nb(2))   1 - (beta/w)*1i*linspace(0,1,nb(2)).^l]';
    
    if nb(3)-2 >= sum(Part(1:id))
        p3 = [1 - (beta/w)*1i*linspace(1,0,nl(3)).^l]';
    else if sum(Part(1:id-1)) < nb(3)+2
            if id == 1
                nbl = nb(3);
            else
                nbl = nb(3)+2-sum(Part(1:id-1));
            end
            p3  = [1 - (beta/w)*1i*linspace(1,0,nbl).^l, ones(1,nl(3)-nbl)]';
        else if nb(3)-2 >=sum(Part(id:end))
                p3 = [1 - (beta/w)*1i*linspace(1,0,nl(3)).^l]';
            else if sum(Part(id+1:end)) < nb(3)+2
                    if id == nlab
                        nbl = nb(3);
                    else
                        nbl = nb(3)+2-sum(Part(id+1:end));
                    end
                    p3  = [ones(1,nl(3)-nbl), 1 - (beta/w)*1i*linspace(0,1,nbl).^l]';
                else
                    p3  = [ones(1,nl(3))]';
                end
            end
        end
    end
    
    %% End of modification
    
    p1 = repmat(permute(p1,[1 2 3]),[1 nl(2) nl(3)]); p1 = p1(:);
    p2 = repmat(permute(p2,[3 1 2]),[nl(1) 1 nl(3)]); p2 = p2(:);
    p3 = repmat(permute(p3,[2 3 1]),[nl(1) nl(2) 1]); p3 = p3(:);
    
    %% averaging matrices
    %%%%%%%%%
    
    %% allocote memory
    R  = zeros(prod(nl),27);
    
    %% Mass matrix
    wm1 = 0.4964958;
    wm2 = 0.4510125;
    wm3 = 0.052487;
    wm4 = 0.45523e-5;
    
    
    R = stencil_add(R,[0 0 0],wm1,w^2*mln.*bln);
    I = [1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];
    R = stencil_add(R,I,(wm2/6),shift(w^2*mln.*bln,I,nl));
    I = [1 1 0;0 1 1;1 0 1;-1 1 0;0 -1 1;-1 0 1;1 -1 0;0 1 -1;1 0 -1;-1 -1 0;0 -1 -1;-1 0 -1];
    R = stencil_add(R,I,(wm3/12),shift(w^2*mln.*bln,I,nl));
    I = [1 1 1;-1 -1 -1;-1 1 1;1 -1 1;1 1 -1;-1 -1 1;1 -1 -1;-1 1 -1];
    R = stencil_add(R,I,(wm4/8),shift(w^2*mln.*bln,I,nl));
    
    %% Stifness matrix
    w1 = 1.8395265e-5;
    w2 = 0.890077;
    w3 = 0.1099046;
    
    
    %%%%%%%%%%%%%%%
    
    % Sc,  eq. (D-2)
    Ix = [1 0 0;-1  0  0];
    Iy = [0 1 0; 0 -1  0];
    Iz = [0 0 1; 0  0 -1];
    
    
    R = getSx(R,Ix,nl,bln,d,p1,w1);
    R = getSy(R,Iy,nl,bln,d,p2,w1);
    R = getSz(R,Iz,nl,bln,d,p3,w1);
    % Sx
    Ix = [ 1  0  0;-1  0  0];
    Iy = [ 0  1  1; 0 -1 -1; 0  1 -1; 0 -1  1];
    Iz = [ 0  1  1; 0 -1 -1; 0 -1  1; 0  1 -1];
    R =  getSx(R,Ix,nl, bln,d,p1,w2/3);
    R =  getSy(R,Iy,nl, bln,d,p2,.25*(w2/3));
    R =  getSz(R,Iz,nl, bln,d,p3,.25*(w2/3));
    % Sy
    Ix = [ 1  0  1;-1  0 -1; 1 0 -1;-1 0  1];
    Iy = [ 0  1  0; 0 -1  0];
    Iz = [ 1  0  1;-1  0 -1;-1 0  1; 1 0 -1];
    R = getSx(R,Ix,nl, bln,d,p1,.25*(w2/3));
    R = getSy(R,Iy,nl, bln,d,p2,w2/3);
    R = getSz(R,Iz,nl, bln,d,p3,.25*(w2/3));
    % Sz
    Ix = [ 1  1  0;-1 -1  0; 1 -1  0;-1  1  0];
    Iy = [ 1  1  0;-1 -1  0;-1  1  0; 1 -1  0];
    Iz = [ 0  0  1; 0  0 -1];
    R = getSx(R,Ix,nl, bln,d,p1,.25*(w2/3));
    R = getSy(R,Iy,nl, bln,d,p2,.25*(w2/3));
    R = getSz(R,Iz,nl, bln,d,p3,w2/3);
    % S1
    Ix = [ 1  1 -1;-1 -1  1; 1 -1  1;-1  1 -1];
    Iy = [ 1  1  1;-1 -1 -1;-1  1 -1; 1 -1  1];
    Iz = [ 1  1  1;-1 -1 -1;-1 -1  1; 1  1 -1];
    R = getSx(R,Ix,nl, bln,d,p1,.25*(w3/4));
    R = getSy(R,Iy,nl, bln,d,p2,.25*(w3/4));
    R = getSz(R,Iz,nl, bln,d,p3,.25*(w3/4));
    % S2
    Ix = [ 1  1  1;-1 -1 -1; 1 -1 -1;-1  1  1];
    Iy = [ 1  1 -1;-1 -1  1;-1  1  1; 1 -1 -1];
    Iz = [ 1  1  1;-1 -1 -1;-1 -1  1; 1  1 -1];
    R = getSx(R,Ix,nl, bln,d,p1,.25*(w3/4));
    R = getSy(R,Iy,nl, bln,d,p2,.25*(w3/4));
    R = getSz(R,Iz,nl, bln,d,p3,.25*(w3/4));
    % S3
    Ix = [ 1  1  1;-1 -1 -1; 1 -1 -1;-1  1  1];
    Iy = [ 1  1  1;-1 -1 -1;-1  1 -1; 1 -1  1];
    Iz = [ 1 -1  1;-1  1 -1;-1  1  1; 1 -1 -1];
    R = getSx(R,Ix,nl, bln,d,p1,.25*(w3/4));
    R = getSy(R,Iy,nl, bln,d,p2,.25*(w3/4));
    R = getSz(R,Iz,nl, bln,d,p3,.25*(w3/4));
    % S4
    Ix = [ 1  1 -1;-1 -1  1; 1 -1  1;-1  1 -1];
    Iy = [ 1  1 -1;-1 -1  1;-1  1  1; 1 -1 -1];
    Iz = [ 1 -1  1;-1  1 -1;-1  1  1; 1 -1 -1];
    R = getSx(R,Ix,nl, bln,d,p1,.25*(w3/4));
    R = getSy(R,Iy,nl, bln,d,p2,.25*(w3/4));
    R = getSz(R,Iz,nl, bln,d,p3,.25*(w3/4));
    %%
    B  = kron(Rf(nl(3)), kron(Rf(nl(2)), Rf(nl(1))));
    e  = B*ones(prod(nl),1);
    I = find(e==0);
    
    idx = [];
    
    % At this line we need the global n
    for k=[-1:1],for l=[-1:1],for mm=[-1:1] idx = [idx mm + nl(1)*l + nl(1)*nl(2)*k];end,end,end
    
    for k=[1:13 15:27]
        Ik = mod(I-1-idx(k),prod(nl))+1;
        R(Ik,k) = 0;
    end
    
    
    if id == 1
        R(end-2*nl(1)*nl(2)+1:end,:) = [];
    else if id == nlab
            R(1:2*nl(1)*nl(2),:) = [];
        else
            R(1:2*nl(1)*nl(2),:) = [];
            R(end-2*nl(1)*nl(2)+1:end,:) =[];
        end
    end
    
    
    for i = 1:length(Part)
        Partition(i) = n(1)*n(2)*Part(i);
    end
    
    codistr = codistributor1d(1, Partition, [prod(n),27]);
    RR       = codistributed.build(R,codistr,'noCommunication');
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Af(k,n)
if k==1
    A = spdiags(ones(n,1)*[.5 .5],[0 1],n,n);
else if k==-1
        A = spdiags(ones(n,1)*[.5 .5],[-1 0],n,n);
    else
        A = spdiags(ones(n,1),0,n,n);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function S = getS(S,b,I,n)
%
% N = prod(n);
% m = size(I,1);
% w = sign(mod(1:m,2)-.5);
%
% for k = 1:m
%     Ik = I/2 + ones(m,1)*I(k,:)/2;
%     ak = w(k)*b(I(k,:));
%     S = stencil_add(S,Ik,ak,w);
% end






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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Av(k,nl)
y = opKron(Af(k(3),nl(3)), Af(k(2),nl(2)),Af(k(1),nl(1)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = bx(k,n,bln,d,p1)
y = (Av(k,n)*bln)./(Av(k,n)*p1)./(d(1).^2*p1);



function y = by(k,n,bln,d,p2)
y = (Av(k,n)*bln)./(Av(k,n)*p2)./(d(2).^2*p2);



function y = bz(k,n,bln,d,p3)
y = (Av(k,n)*bln)./(Av(k,n)*p3)./(d(3).^2*p3);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Rf(n)
y = spdiags([0; ones(n-2,1); 0],0,n,n);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = getSx(S,I,n,bln,d,p,ww)

N = prod(n);
m = size(I,1);
w = sign(mod(1:m,2)-.5);

for k = 1:m
    Ik = I/2 + ones(m,1)*I(k,:)/2;
    ak = w(k)*ww*bx(I(k,:),n,bln,d,p);
    S  = stencil_add(S,Ik,ak,w);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = getSy(S,I,n,bln,d,p,ww)

N = prod(n);
m = size(I,1);
w = sign(mod(1:m,2)-.5);

for k = 1:m
    Ik = I/2 + ones(m,1)*I(k,:)/2;
    ak = w(k)*ww*by(I(k,:),n,bln,d,p);
    S  = stencil_add(S,Ik,ak,w);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = getSz(S,I,n,bln,d,p,ww)

N = prod(n);
m = size(I,1);
w = sign(mod(1:m,2)-.5);

for k = 1:m
    Ik = I/2 + ones(m,1)*I(k,:)/2;
    ak = w(k)*ww*bz(I(k,:),n,bln,d,p);
    S  = stencil_add(S,Ik,ak,w);
end























