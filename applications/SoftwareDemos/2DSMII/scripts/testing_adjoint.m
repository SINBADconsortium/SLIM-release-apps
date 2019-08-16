% This script does a `dot-test' of the Jacobian for a small model.
% It should take only a couple of minutes to run.s
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
model.n = [51 51];
model.nb = [10 10];
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

%
J = oppDF(m,Q,model);

% test for 10 random vectors
a = zeros(10,1);
b = zeros(10,1);
for k = 1:10
    x = randn(prod(model.n),1);
    y = J*randn(prod(model.n),1); 
    
    a(k) = real(gather((J*x)'*y));
    b(k) = real(gather(x'*(J'*y)));
end

% save
dlmwrite('error_adj.dat',[a(:) b(:) abs(a(:)-b(:)).^2]);

cd(curdir)