% This script tests the modeling operator against the analytic solution
% for a constant velocity. It should take a couple of minutes to run.

setpaths;
curdir = pwd;
expdir = [resultsdir '/testing'];
if ~exist(expdir,'dir')
    mkdir(expdir);
end
cd(expdir);

% Analytic solution for constant velocity
% velocity [m/s]
v0 = 2000;

% size of domain LxL m
L = 1000;

% sponge boundary size B m
B = 250;

% frequency [1/s]
f  = 10;

% gridspacings
h = [20:-2:2];
e = 0*h;
% loop
for k = 1:length(h)
    % grid paramets
    model.o = [0 0];
    model.d = [h(k) h(k)];
    model.n = [floor(L/h(k)) + 1 floor(L/h(k)) + 1];

    % source and receiver setup
    model.zsrc = 500;
    model.xsrc = 500;
    model.zrec = 0:h(k):L;
    model.xrec = 0:h(k):L;

    % absorbing boundary
    model.nb = [floor(B/h(k)) + 1 floor(B/h(k))];

    % frequency
    model.freq = f;

    % ricker wavelet
    model.f0 = 10;
    model.t0 = 0;

    % source matrix
    Q = 1;

    % slowness-squares model [km^2/s^2]
    m  = 1e6*ones(prod(model.n),1)./v0.^2;

    % analytic solution
    D1 = G([v0;0],Q,model);

    % numerical solution
    D2 = F(m,Q,model);

    % error
    e(k) = norm(D1 - D2);
end

% save
dlmwrite('error_mod.dat',[h(:) e(:)]);
[o,d,n] = grid2odn(model.zrec,model.xrec);
odnwrite('D1.odn',D1(:),o,d,n);
odnwrite('D2.odn',gather(D2(:)),o,d,n);


% Analytic solution for linear profile
% velocity profile
v0 = 2000;
alpha = 0.7;

% grid
z = 0:10:1000;
x = 0:10:5000;
[o,d,n] = grid2odn(z,x);

% gridded slowness-squared
m = vec((1e6./(v0+alpha*z').^2)*ones(size(x)));

% frequency
f = 10;

% modeling parameters
model.o = o;
model.d = d;
model.n = n;
model.nb = [50 50];
model.freq = f;
model.f0 = 10;
model.t0 = 0;
model.zsrc = 10;
model.xsrc = mean(x);
model.zrec = z;
model.xrec = x;

% source
Q = 1;

% analytic solution
G1 = G([v0;alpha],Q,model);

% numerical solution
G2 = F(m,Q,model);

% save
[o,d,n] = grid2odn(model.zrec,model.xrec);
odnwrite('G1.odn',G1(:),o,d,n);
odnwrite('G2.odn',gather(G2(:)),o,d,n);

cd(curdir)
