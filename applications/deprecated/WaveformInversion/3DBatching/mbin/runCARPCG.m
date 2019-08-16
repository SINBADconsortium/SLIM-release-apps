function runCARPCG(expdir,vfile,p,ns,bs,tol,maxit)
% run CARPCG
% 
% use:
%   runCARPCG(vfile,p,ns,bs,tol,maxit)
%
% input:
%   vfile - velocity model in .odn format
%   p     - subampling in each direction
%   ns    - number of sources
%   bs    - number of blocks (for block-cg)
%   tol   - tolerance
%   maxit - max. number of iterations

if ~exist(expdir,'dir')
    mkdir(expdir);
end
curdir = pwd;
cd(expdir);

spmd, np = numlabs; end;
np = np{1};

if np>1
    ns = 1;
end

[v,n,d,o] = rsf_read_all(vfile);
vmin = min(v(:));
v  = reshape(v,n);
v  = v(1:p:end,1:p:end,1:p:end);
d  = d*p;
f  = .5*round(2*vmin/(10*max(d)));
nb = ceil(500./d);

n = size(v);
N = prod(n);
rsf_write_all(['v_' num2str(p) '_' num2str(np) '.rsf'],{'out=stdout'},v,d,o);

% construct source and matrix
beta = 100;
l    = 2;

% get sources
Q  = zeros(N,ns);
Q(sub2ind(n,nb(1) + ceil(50/d(1)) + 1,floor(n(2)/2)+1,floor(n(3)/2)+1)) = -1/prod(d);
for k = 2:ns 
    ix = nb(2) + randi(n(2)-2*nb(2),1,1);
    iy = nb(3) + randi(n(3)-2*nb(3),1,1);
    Q(sub2ind(n,nb(1) + ceil(50/d(1)) + 1,ix,iy),k) = -1/prod(d);
end

[R,idx] = R_Helm3D(f,1./v(:).^2,ones(N,1),o,d,n,nb,beta,l);

% normalize matrix
a = 1./sqrt(sum(abs(R).^2,2));
R = bsxfun(@times,R,a);
Q = bsxfun(@times,Q,a);

clear v a;

% distribute matrix
if np > 1
    x0 = distributed.zeros(n); x0 = x0(:);
    Q  = distributed(Q);
    R  = distributed(R);
    spmd,
        dist1 = getCodistributor(x0);
        dist2 = codistributor1d(1,dist1.Partition,[prod(n),27]);
        Q    = redistribute(Q,dist1);
        R    = redistribute(R,dist2);
    end
else
    x0 = zeros(N,ns);
end

% time CG
options.tol   = tol;
options.maxit = maxit;
options.w     = 1.5;
options.nb    = bs;
for k = 1:1
    tic
    if np==1
        % This actually runs "block" CGMN, and not CARPCG
        [u,res,niter] = CARPBCG(R,idx,Q,x0,n,options);
    else
        [u,res] = CGMN(R,idx,Q,x0,options);
        niter = length(res);
    end
    Tcg(k) = toc;
end
Tcg = mean(Tcg);

% output
if exist(['table.dat'],'file')
    fid = fopen(['table.dat'],'a');
else
    fid = fopen(['table.dat'],'w');
    fprintf(fid,'f, ng, dx, N, np, ns, bs, niter, Tcg\n');
end   
fprintf(fid,'%2.1f, %2.1f, %2d, %6d, %2d, %2d, %2d, %5d, %4.1f\n',f, vmin/(f*max(d)),max(d), prod(n),np,ns,bs,niter,Tcg);
fclose(fid);

if np>1
    u = gather(u);
end

% write results
rsf_write_all(['u_' num2str(p) '_' num2str(np) '.rsf'],{'out=stdout'},reshape(u(:,1),n),d,o);
dlmwrite(['res_' num2str(p) '_' num2str(np) '_' num2str(bs) '.dat'],res);

cd(curdir);
