function test_suite = testWRI
% unit tests for Waveform Reconstruction Inversion
%
% Curt Da Silva, 2015
% curtd@math.ubc.ca
%
initTestSuite;

function model = setup
model.o = [0 0 0];
model.d = [10 10 1]; 
nx = 151; nz = 151;
model.n = [nz nx 1]; 
model.nb = [40 40 0];
if parpool_size() > 0
    model.freq = linspace(10,15,parpool_size());
else
    model.freq = [10 15];
end
model.unit = 's2/km2';
model.f0 = 10;
model.t0 = 0;
model.zsrc = 50;
model.xsrc = 100:100:1000;
model.zrec = 50;
model.xrec = 0:10:1500;
x2 = linspace(0,1,nx); z2 = linspace(0,1,nz);
[Z,X] = ndgrid(z2,x2);
m = 1e6/2000^2 * ones(nz,nx);
m( 8*(Z-0.4) + X >= 0 ) = 1e6/2500^2;
m( 8*(Z-0.8) + X >= 0 ) = 1e6/3000^2;
model.m = m;


function testGradient(model)
% Ensures that the WRI gradient passes the gradient test

% Rate of convergence tolerance (due to inversion of matrices, floating point errors may inhibit true convergence when the step size is small)
tol = 0.2;

% physical grid
[z,x]  = odn2grid(model.o,model.d,model.n);
m = model.m;
% comp. grid
dm = randn(model.n); dm(1,:) = 0; dm(end,:) = 0; dm(:,1) = 0; dm(:,end) = 0;
dm = dm/norm(vec(dm))*1e-1*max(abs(vec(m)));
m = m(:); dm = dm(:);
nsrc = length(model.xsrc);
Q = speye(nsrc);

Dobs = F(m+dm,Q,model);
params = struct;
params.lambda = 1e1;
params = default_fwi_params2d(params);
[f0,g0] = misfit_pen(m,Q,Dobs,model,params);
if parpool_size()==0
    H = opHpen(m,Q,Dobs,model,params);
else
    H = oppHpen(m,Q,Dobs,model,params);
end
dm = randn(model.n); dm(1,:) = 0; dm(end,:) = 0; dm(:,1) = 0; dm(:,end) = 0;
dm = dm/norm(vec(dm))*1e-1*max(abs(vec(m))); dm = dm(:);

h0dm = H*dm;
obj = @(h) misfit_pen(m+h*dm,Q,Dobs,model,params);
dobj = g0' * dm;
d2obj = dm'*h0dm;

h = 10.^(-3:2);
e0 = zeros(length(h),1);
e1 = zeros(length(h),1);
e2 = zeros(length(h),1);
for i=1:length(h)
    fh = obj(h(i));
    e0(i) = abs(fh - f0);
    e1(i) = abs(fh - f0 - h(i)*dobj);
    e2(i) = abs(fh - f0 - h(i)*dobj - 0.5*h(i)^2*d2obj);
end
h0 = median(log10(e0(2:end)./e0(1:end-1)));
h1 = median(log10(e1(2:end)./e1(1:end-1)));
h2 = median(log10(e2(2:end)./e2(1:end-1)));

% First order convergence without gradient
assert( abs(h0 - 1) < tol);

% Second order convergence with gradient
assert( abs(h1 - 2) < tol);

% Third order convergence with hessian
%assert( abs(h2 - 3) < tol);

%Dot test
y = randn(size(H,2),1); y = reshape(y,nz,nx); y([1 end],:) = 0; y(:,[1 end]) = 0; y = vec(y);
z = randn(size(H,2),1); z = reshape(z,nz,nx); z([1 end],:) = 0; z(:,[1 end]) = 0; z = vec(z);

s = (H*y); s = z'*s;
t = (H*z); t = y'*t;

assert(abs(s - t) < 1e-10 * max(abs(s),abs(t)));
