% This script runs a Jacobian test on a small model.
% It should take only a couple of minutes to run.
%

setpaths;
curdir = pwd;
expdir = [resultsdir '/testing'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% setup model parameters
model.o = [0 0];
model.d = [10 10];
model.n = [101 101];
model.nb = [20 20];
model.freq = [10 15];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 15;
model.xsrc = 0:100:1000;
model.zrec = 10;
model.xrec = 0:5:1000;

% source matrix
Q = speye(length(model.xsrc));

% constant velocicty
v0 = 2000;
m  = 1e6/v0.^2*ones(prod(model.n),1);

% random perturbation, set values close to the edge to zero.
dm = randn(model.n);
dm([1:20 end-20:end],:) = 0;
dm(:,[1:20 end-20:end]) = 0;
dm = dm(:);
dm = dm/max(dm)*1e-1*max(m);

% data and Jacobian
[D0,J0] = F(m,Q,model);

% linearized data
dD = J0*dm;

% stepsizes
h = 10.^[0:-1:-6];
e = 0*h;
for k = 1:length(h)
    D1   = F(m+h(k)*dm,Q,model);
    e(k) = norm(D1 - D0 - h(k)*dD)/norm(dm);
end

% save
dlmwrite('error_jac.dat',[h(:) e(:)]);


cd(curdir)