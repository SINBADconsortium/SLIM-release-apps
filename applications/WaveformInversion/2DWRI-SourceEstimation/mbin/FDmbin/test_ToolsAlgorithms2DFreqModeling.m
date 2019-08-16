function test_suite = test_ToolsAlgorithms2DFreqModeling
% unit tests for 2DFDFDModeling
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% Updated by
% Curt Da Silva, 2015
% curtd@math.ubc.ca
%
initTestSuite;

function model = setup
model.o = [0 0 0];
model.d = [10 10 1];
model.n = [101 101 1];
model.nb = [20 20 0];
if parpool_size() > 1
    model.freq = linspace(4,8,parpool_size());
else
    model.freq = 4;
end
model.f0 = 10;
model.t0 = 0.0;
model.zsrc = 15;
model.xsrc = 500;
model.zrec = 10;
model.xrec = 0:5:1000;
model.unit = 's2/km2';
v0 = 2000;
model.m = 1e6/v0.^2*ones(prod(model.n),1);
rng(12385352);

function testHelmholtz(model)
m = model.m;
model.unit = 's2/km2';

dt = model.d;
nt = model.n+2*model.nb;
Nt = prod(nt);
model.nb = repmat(model.nb,2,1);
dm = randn(model.n);

Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
tol = 0.2;
% scan over frequencies
freq = 1:10;
m = m(:); dm = dm(:);
mx = Px*m; dmx = Px*dm;
h = 10.^(-6:0);

for i=1:length(freq)
    [H0,dH0,ddH0] = Helm2D_opt(mx,dt,nt,model.nb(:,1:2),model.unit,freq(i),model.f0);
    e0 = zeros(size(h)); e1 = zeros(size(h)); e2 = zeros(size(h));
    dH0 = dH0 * spdiags(dmx,0,Nt,Nt);
    ddH0 = ddH0 * spdiags(dmx.^2,0,Nt,Nt);
    for j=1:length(h)
        H1 = Helm2D_opt(mx+h(j)*dmx,dt,nt,model.nb(:,1:2),model.unit,freq(i),model.f0);
        e0(j) = norm(H1-H0,'fro');
        e1(j) = norm(H1-H0-h(j)*dH0,'fro');
        e2(j) = norm(H1-H0-h(j)*dH0 - h(j)^2/2 * ddH0,'fro');
    end
    h0 = log10(mean(e0(2:end)./e0(1:end-1)));
    assert( abs( h0 - 1) < tol );
    assert( max(e1) < 1e-10 );
    assert( max(e2) < 1e-10 );
end

model.unit = 'm/s'; mx = 1e3*(mx.^(-1/2));
for i=1:length(freq)
    [H0,dH0,ddH0] = Helm2D_opt(mx,model.d,nt,model.nb(:,1:2),model.unit,freq(i),model.f0);
    e0 = zeros(size(h)); e1 = zeros(size(h)); e2 = zeros(size(h));
    dH0 = dH0 * spdiags(dmx,0,Nt,Nt);
    ddH0 = ddH0 * spdiags(dmx.^2,0,Nt,Nt);
    for j=1:length(h)
        H1 = Helm2D_opt(mx+h(j)*dmx,model.d,nt,model.nb(:,1:2),model.unit,freq(i),model.f0);
        e0(j) = norm(H1-H0,'fro');
        e1(j) = norm(H1-H0-h(j)*dH0,'fro');
        e2(j) = norm(H1-H0-h(j)*dH0 - h(j)^2/2 * ddH0,'fro');
    end
    h0 = log10(median(e0(2:end)./e0(1:end-1)));
    h1 = log10(median(e1(2:end)./e1(1:end-1)));
    
    assert( abs( h0 - 1) < tol );
    assert( abs( h1 - 2) < tol );
    assert( max(e2) < 1e-8 );
end
%m = 1e3*(m.^(-1/2));
model.unit = 's2/km2';
params = struct;
params.pml = 20;
params.mat_free = false;
params.scheme = PDEopts2D.default_helm_scheme;
params.solve_opts = LinSolveOpts2D();
[H,~,dH,ddH] = discrete_helmholtz(m,model,freq(1),params);
dH = dH * spdiags(dmx,0,Nt,Nt);
ddH = ddH * spdiags(dmx.^2,0,Nt,Nt);
e0 = zeros(size(h)); e1 = zeros(size(h)); e2 = zeros(size(h));
for j=1:length(h)
    H1 = discrete_helmholtz(m+h(j)*dm,model,freq(1),params);
    e0(j) = norm(H1.params.coef-H.params.coef,'fro');
    e1(j) = norm(H1.params.coef-H.params.coef-h(j)*dH,'fro');
    e2(j) = norm(H1.params.coef-H.params.coef-h(j)*dH - h(j)^2/2 * ddH,'fro');
end
h0 = log10(median(e0(2:end)./e0(1:end-1)));

assert( abs( h0 - 1) < 0.2 );
assert( abs( h1 - 2) < 0.2 );
assert( max(e2) < 1e-8 );

function testAnalytic(model)
model.n(1:2) = 2*model.n(1:2);
model.d(1:2) = model.d(1:2)/2;
nsrc = length(model.xsrc);
Q = speye(nsrc);
v0 = 2000;
vc = v0*ones(model.n);
freq = model.freq(1);
model.unit = 'm/s';

Ptr = opKron(opExtension(model.n(2),[model.nb(2) model.nb(2)],0),opExtension(model.n(1),[model.nb(1) model.nb(1)],0));
params = struct;
params.pml = 20;
params.mat_free = false;
params.scheme = PDEopts2D.default_helm_scheme;
params.solve_opts = LinSolveOpts2D();
[Hk,comp_grid] = discrete_helmholtz(vc,model,freq,params);
[zt,xt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc))';
q = Ps*Q(:,1);
u = Hk\q;
[~,I] = max(abs(vec(q))); [iz,ix] = ind2sub(comp_grid.nt,I);
[ZT,XT] = ndgrid(zt,xt);
R = ((ZT-zt(iz)).^2 + ((XT-xt(ix)).^2)).^(1/2);
t = (2*pi*freq/v0)*R;
G = -1i/4*besselh(0,t);
G(iz,ix) = -0.5+0.25*1i;
G = prod(model.d)* conj(G);

r = G-reshape(u,comp_grid.nt);
r(iz,ix) = 0; %We don't expect an accurate wavefield at the source
assert(norm(Ptr'*vec(r))/norm(Ptr'*vec(G)) < 0.01);


function testdU(model)
%Test derivative of forward wavefield
m = model.m;
% physical grid
[z,x]  = odn2grid(model.o,model.d,model.n);

% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
Nt = prod(nt);
[zt,xt] = odn2grid(ot,dt,nt);
nz = model.n(1); nx = model.n(2);
dm = randn(nz,nx);
model.nb = repmat(model.nb,2,1);

Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));

m = Px*m(:); dm = Px*dm(:);

% Ricker wavelet
w = fwi_wavelet(model.freq,model.t0,model.f0);
Q = speye(length(model.xsrc));
h = 10.^(-6:0);

tol = 0.15;
i=1;
[H0,dH0] = Helm2D_opt(m,dt,nt,model.nb(:,1:2),model.unit,model.freq(i),model.f0);
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
q = Ps'*(w(i) * Q);
U0 = H0\q;
dHk = dH0*spdiags(dm,0,Nt,Nt);
dU = H0\(dHk*(-U0));
e0 = zeros(length(h),1);
e1 = zeros(length(h),1);
for k=1:length(h)
    m1 = m + h(k)*dm;
    H = Helm2D_opt(m1,dt,nt,model.nb(:,1:2),model.unit,model.freq(i),model.f0);
    U1 = H\q;
    e0(k) = norm(U1-U0,'fro');
    e1(k) = norm(U1-U0-h(k)*dU,'fro');
end
% order of convergence
h0 = log10(mean(e0(2:end)./e0(1:end-1)));
h1 = log10(mean(e1(2:end)./e1(1:end-1)));

% first order convergence without gradient
assert( abs(h0 - 1 ) < tol );
% second order convergence with gradient
assert( abs(h1 - 2 ) < tol );

% Test adjoint relationship
Z = randn(size(U0))+1i*randn(size(U0));
x = vec(dU)' * vec(Z);
V = H0'\(-Z);
y = vec(dm)' * vec( sum( (conj(U0) .* (dH0'*V)) , 2) );
assert(abs(x-y) < 1e-10*max(abs(x),abs(y)));


function testLU(model)
% Ensures that forward modeling, jacobian + transpose
% produce the same answer when the LU decomposition flag is turned
% on or off
m = model.m;
Q = speye(length(model.xsrc));
solve_opts = LinSolveOpts;
solve_opts.solver = LinSolveOpts.SOLVE_BACKSLASH;
params.lsopts = solve_opts;
[D1,J1] = F(m,Q,model,params);

solve_opts.solver = LinSolveOpts.SOLVE_LU;
params.lsopts = solve_opts;

[D2,J2] = F(m,Q,model,params);

assertVectorsAlmostEqual(gather(D1),gather(D2),'relative');

dm = randn(size(m));

J1dm = J1*dm; J2dm = J2*dm;

assertVectorsAlmostEqual(gather(J1dm),gather(J2dm),'relative');

J1t = J1'*J1dm; J2t = J2'* J2dm;

assertVectorsAlmostEqual(gather(J1t),gather(J2t),'relative');


function testSource(model)
% test frequency-dependent source
% added: Feb. 1, 2012
model.xsrc = 0:10:1000;
model.xrec = 0:10:1000;
nsrc = length(model.xsrc);
nrec = length(model.xrec);
nfreq = length(model.freq);

Q = speye(nsrc);
m = 1e6*ones(prod(model.n),1)./2000.^2;

D  = F(m,Q,model);
D  = invvec(D,[nrec nsrc nfreq]);
DD = F(m,D,model);
DD = invvec(DD,[nrec nsrc nfreq]);

D = gather(D); DD = gather(DD);
for k = 1:nfreq
    tmpa = vec(D(:,:,k)*D(:,:,k));
    tmpb = vec(DD(:,:,k));
    assertVectorsAlmostEqual(tmpa,tmpb,'relative',0.05);
end



function testJacobian1(model)
% test forward Jacobian taylor error
Q = speye(length(model.xsrc));
m = model.m;
[D0, J0] = F(m,Q,model);
D0 = gather(D0);

%Randomize over perturbations
nTrials = 1;

%Error in order of convergence tolerance, since the order of convergence is asymptotic + limited precision from having to solve linear systems
tol = 0.15;

% if h << 1e-6, you run in to floating point issues so the test might fail then
h = 10.^(-6:0);

for j=1:nTrials
    
    dm = randn(model.n); dm(1,:) = 0; dm(end,:) = 0;
    dm(:,1) = 0; dm(:,end) = 0;
    dm = dm(:); dm = dm/max(dm)*1e-1*max(m);
    
    dD = gather( J0*dm );
    
    % Gradient test
    e0 = zeros(length(h),1);
    e1 = zeros(length(h),1);
    for i=1:length(h)
        D1 = gather(F(m+h(i)*dm,Q,model));
        e0(i) = norm(D1-D0,'fro');
        e1(i) = norm(D1-D0-h(i)*dD,'fro');
    end
    
    % Order of convergence
    h0 = log10(mean(e0(2:end)./e0(1:end-1)));
    h1 = log10(mean(e1(2:end)./e1(1:end-1)));
    
    % First order convergence without gradient
    assert( abs( h0 - 1 ) < tol );
    
    % Second order convergence with gradient
    assert( abs( h1 - 2 ) < tol );
end


function testJacobian2(model)
% test Jacobian with extended source
% added: Feb 1, 2012
Q = speye(length(model.xsrc));
m = model.m;

tol = 0.15;

dm = randn(model.n);
dm = dm(:); dm = dm/max(dm)*1e-1*max(m);
if parpool_size()>0
    Q = distributed.randn(length(model.xsrc),length(model.xsrc),length(model.freq));
else
    Q = randn(length(model.xsrc),length(model.xsrc),length(model.freq));
end
[D0, J0] = F(m,Q,model);
dD = J0*dm;

h = 10.^(-6:0);
e0 = zeros(length(h),1);
e1 = zeros(length(h),1);
for i=1:length(h)
    D1 = gather(F(m+h(i)*dm,Q,model));
    e0(i) = norm(D1-D0,'fro');
    e1(i) = norm(D1-D0-h(i)*dD,'fro');
end

% Order of convergence
h0 = log10(mean(e0(2:end)./e0(1:end-1)));
h1 = log10(mean(e1(2:end)./e1(1:end-1)));

% First order convergence without gradient
assert( abs( h0 - 1 ) < tol );

% Second order convergence with gradient
assert( abs( h1 - 2 ) < tol );


function testJacobian3(model)
% Test LS objective gradient convergence
Q = speye(length(model.xsrc));
m = model.m;

tol = 0.2;

dm = randn(model.n);
dm([1 end],:) = 0; dm(:,[1 end]) = 0;
dm = dm(:); dm = dm/max(dm)*1e-1*max(m);

[D0, J0] = F(m,Q,model);
D1 = F(m+dm,Q,model);
obj = @(h) 0.5* norm(F(m+h*dm,Q,model)-D1,'fro')^2;
r = D0-D1; f0 = 0.5 * norm(r,'fro')^2;
g0 = J0' * r;
dobj = g0' * (dm);
h = 10.^(-6:0);
e0 = zeros(length(h),1);
e1 = zeros(length(h),1);

for i=1:length(h)
    fh = obj(h(i));
    e0(i) = abs(fh-f0);
    e1(i) = abs(fh - f0 - h(i)*dobj);
end
h0 = log10(mean(e0(2:end)./e0(1:end-1)));
h1 = log10(mean(e1(2:end)./e1(1:end-1)));

% First order convergence without gradient
assert(h0 > 1 || 1-h0 < tol );

% Second order convergence with gradient
assert(abs(h1 - 2) < tol );

function testAdjoint1(model)
% Adjoint test of Jacobian
Q = speye(length(model.xsrc));
m = model.m;

dm = randn(model.n); dm([1 end],:) = 0; dm(:,[1 end]) = 0; dm = vec(dm);
dm = dm/max(dm)*1e-1*max(m);

Q = speye(length(model.xsrc));

[~,J0] = F(m,Q,model);
dD = J0*dm;
Z = randn(size(dD));
if isdistributed(dD)
    Z = distributed(Z);
    spmd, Z = redistribute(Z,getCodistributor(dD)); end
end

A = real(gather(vec(dD)' * vec(Z)));
B = vec(dm)' * vec(J0'*Z);

assert( abs(A - B) < 1e-10*max(abs(A),abs(B)) );


function testAdjoint2(model)
% Adjoint test of Jacobian with extended source
Q = speye(length(model.xsrc));
m = model.m;

dm = randn(model.n); dm([1 end],:) = 0; dm(:,[1 end]) = 0; dm =vec(dm);
dm = dm/max(dm)*1e-1*max(m);
Q  = distributed(randn(length(model.xsrc),length(model.xsrc),length(model.freq)));

J0 = oppDF(m,Q,model,[]);
dD = J0*dm;
Z = randn(size(dD));
if isdistributed(dD)
    Z = distributed(Z);
    spmd, Z = redistribute(Z,getCodistributor(dD)); end
end

A  = real(gather(dD'*Z));
B  = dm'*(J0'*Z);

assertElementsAlmostEqual(gather(A),gather(B));


function testH(model)
nsrc = length(model.xsrc);
nrec = length(model.xrec); nfreq = length(model.freq);
freq = model.freq;
Q = speye(nsrc);
m = model.m;
opts = struct;
opts = default_fwi_params2d(opts);
dm = randn(model.n);
dm([1 end],:) = 0; dm(:,[1 end]) = 0;
dm = vec(dm); dm = dm/max(dm)*1e-1*max(m);
Dobs = F(m+dm,Q,model); 
dm = randn(model.n);
dm(1,:) = 0; dm(end,:) = 0; dm(:,1) = 0; dm(:,end) = 0;
dm = vec(dm); dm = dm/max(dm)*1e-1*max(m);
if parpool_size()==0
    obj = @(h) PDEfunc(PDEopts.OBJ,m+h*dm,Q,[],Dobs,model,opts);
else
    obj = @(h) PDEfunc_dist(PDEopts.OBJ,m+h*dm,Q,[],Dobs,model,opts);
end
[f0,g0] = obj(0);
Dobs = reshape(Dobs,nrec,nsrc,nfreq);
tol = 0.3;
% Brute force Taylor error calculations
k=1;
D1 = gather(Dobs(:,1,k));
w = fwi_wavelet(freq(k),model.t0,model.f0);
params = struct;
params.pml = 20;
params.mat_free = false;
params.scheme = PDEopts2D.default_helm_scheme;
params.solve_opts = LinSolveOpts2D();
[Hk,comp_grid,dH,ddH] = discrete_helmholtz(reshape(m,model.n),model,freq(k),params);
dHdm = dH*opDiag(comp_grid.phys2comp*dm);
[zt,xt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
Ps = opInterp('sinc',model.zsrc,zt,model.xsrc,xt);
Pr = opInterp('sinc',model.zrec,zt,model.xrec,xt)';

Qk = w(k) * Ps*Q(:,1);

Uk = Hk\Qk;
dUk = Hk\(-dHdm*Uk);
Vk = Hk'\( -Pr'*(Pr*Uk-D1) );
dVk = Hk'\( - dHdm'*Vk - Pr'*Pr*dUk );
g = comp_grid.comp2phys*real(conj(Uk).*(dH'*Vk));
dg = comp_grid.comp2phys*real(conj(dUk) .* (dH'*Vk) + conj(Uk) .* ( ddH'*Vk + dH'*dVk));

h = 10.^(-3:0);
e0 = zeros(size(h)); e1 = zeros(size(h));
for j=1:length(h)
    [H1,~,dH1] = discrete_helmholtz(reshape(m+h(j)*dm,model.n),model,freq(k),params);
    U1 = H1\Qk;
    V1 = H1'\(-Pr'*(Pr*U1-D1));
    g1 = comp_grid.comp2phys*real(conj(U1).*(dH1'*V1));
    
    e0(j) = norm(g1-g);
    e1(j) = norm(g1-g-h(j)*dg);
end
h0 = median(diff(log10(e0)));
h1 = median(diff(log10(e1)));
assert( abs(h0 - 1) < tol || h0 > 1 );
assert( abs(h1 - 2) < tol || h1 > 2 );


% m : s^2/m^2
dm = dm(:);
Dobs = reshape(Dobs,[nsrc*nrec,nfreq]);

if parpool_size()==0
    H = opH(model.m,Q,Dobs,model);
else
    H = oppH(model.m,Q,Dobs,model);
end
Hdm = H*dm;
if parpool_size()==0
    assert(norm(dg-Hdm)<1e-8*norm(Hdm));
end

df = g0'*(dm);
d2f = dm'*Hdm;
h = 10.^(-3:0);
e0 = zeros(size(h));
e1 = zeros(size(h));
e2 = zeros(size(h));
for i=1:length(h)
    f1 = obj(h(i));
    e0(i) = abs(f1-f0);
    e1(i) = abs(f1-f0-h(i)*df);
    e2(i) = abs(f1-f0-h(i)*df-0.5*h(i)^2*d2f);
end

h0 = median(diff(log10(e0)));
h1 = median(diff(log10(e1)));
h2 = max(diff(log10(e2)));

%First, second, third order convergence
assert( abs(h0 - 1) < tol || h0 > 1 );
assert( abs(h1 - 2) < tol || h1 > 2 );
assert( abs(h2 - 3) < tol || h2 > 3 );


%Dot test
y = randn(size(H,2),1); y = reshape(y,model.n); y([1 end],:) = 0; y(:,[1 end]) = 0; y = vec(y);
z = randn(size(H,2),1); z = reshape(z,model.n); z([1 end],:) = 0; z(:,[1 end]) = 0; z = vec(z);
s = (H*y); s = z'*s;
t = (H*z); t = y'*t;
assert(abs(s - t) < 1e-10 * max(abs(s),abs(t)));

