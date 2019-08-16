% 2.5D example of PDEfunc, comparing to the analytic Green's function in a
% constant medium
% 
% Curt Da Silva, 2016

% Number of points in each direction
n1d = 100;

% Number of receivers
nrec = 100;

% Number of sources
nsrc = 100;

% Domain is 2km x 2km x 2km
L = 2000;

% Frequency [hz]
freq = 8;

% Constant velocity [m/s]
v0 = 2000;

% Number of y-wavenumbers to use (i.e., number of 2D PDEs to solve)
% Higher = more accurate, slower
nky = 100;

%% Construct model struct + source weights
model.o = [0 0 0];
model.n = [n1d n1d n1d];
model.d = L./ model.n;
dx = model.d(1); dy = model.d(2); dz = model.d(3);
[x,y,z] = odn2grid(model.o,model.d,model.n);
model.freq = freq;
model.f0 = 10;
model.t0 = 0;
model.unit = 'm/s';
model.zsrc = 20;
model.xsrc = linspace(0,max(x),nsrc);
model.ysrc = 100;
model.zrec = 40;
model.xrec = linspace(0,max(x),nrec);
model.yrec = 0;
nrx = length(model.xrec); nry = length(model.yrec);
Q = speye(length(model.xsrc)*length(model.ysrc));

%% Velocity + green's function
vc = v0*ones(model.n);

lsopts = LinSolveOpts();
lsopts.tol = 1e-10;
lsopts.maxit = 2000;
lsopts.maxinnerit = 5;
lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
lsopts.precond = LinSolveOpts.PREC_MLGMRES;

opts = struct;
opts.free_surface = false;
opts.disp_output = true;
opts.scheme = PDEopts.HELM3D_OPERTO27;
opts.solve_opts = lsopts;
opts.pml = 30;

[Hk,comp_grid] = discrete_helmholtz(vc,model,freq,opts);

[xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);

%% Compare synthesized vs analytic data
vc2D = v0*ones(model.n([3 1]));

opts.pdefunopts = PDEopts2D();
opts.pdefunopts.helm_pml = 30;

% Integration interval is slightly over the critical frequency
%    kc = freq/v0
%
% More accurate for a given number of points, but will cause the Taylor error 
% test to fail because the upper bound depends on the input velocity
opts.use_knyq = false;
opts.nky = nky;
tic,D = PDEfunc_25D(PDEopts.FORW_MODEL, vc2D, Q, [], [],model,opts);
disp(['Data generation time ' num2str(toc) 's']);

% Analytic data 
Dtrue = zeros(length(model.xrec),length(model.xsrc));
Ps = opInterp('sinc',model.xsrc,xt,model.ysrc,yt,model.zsrc,zt);
srcw = fwi_wavelet(model.freq,model.t0,model.f0);
for i=1:nsrc
    q = Ps*Q(:,i);
    [~,I] = max(abs(vec(q))); 
    [ix,iy,iz] = ind2sub(comp_grid.nt,I);
    [XT,YT,ZT] = ndgrid(xt,yt,zt);
    R = ((XT-xt(ix)).^2 +(YT-yt(iy)).^2 + (ZT-zt(iz)).^2 ).^(1/2);
    w = (2*pi*freq/v0);
    G = srcw*prod(model.d)*exp(1i*w*R)./(4*pi*R);
    Gref = squeeze(G(:,yt==model.yrec,:)).';
    Pr = opInterp('sinc',model.zrec,zt,model.xrec,xt)';    
    Dtrue(:,i) = Pr*vec(Gref);
end

isrc = 50;
figure; gp= plot(x,real(Dtrue(:,isrc)),'r','LineWidth',2); hold all;
up=plot(x,real(D(:,isrc)),'b','LineWidth',2); 
legend([gp,up],{'Greens function','2.5D wavefield'});
xlabel('x [m]'); 
title('Real part, 2.5D forward modelled data');

figure; gp= plot(x,imag(Dtrue(:,isrc)),'r','LineWidth',2); hold all;
up=plot(x,imag(D(:,isrc)),'b','LineWidth',2); 
legend([gp,up],{'Greens function','2.5D wavefield'});
xlabel('x [m]'); 
title('Imaginary part, 2.5D forward modelled data');

%% Compute 2.5D Taylor error + adjoint tests

% Reduce number of ky points just for computational purposes
opts.nky = 10;
opts.use_knyq = true;
dv = 100*randn(model_2d.n);
dv([1 end],:) = 0; dv(:,[1 end]) = 0;

tic,D = PDEfunc_25D(PDEopts.FORW_MODEL, vc, Q, [], [],model,opts);disp(toc);
tic,dD = PDEfunc_25D(PDEopts.JACOB_FORW, vc, Q, dv, [],model,opts);disp(toc);

h = 10.^(-6:0);
e0 = zeros(length(h),1); e1 = zeros(length(h),1);
for i=1:length(h)
    D1 = PDEfunc_25D(PDEopts.FORW_MODEL, vc + h(i)*dv, Q, [], [],model,opts);
    e0(i) = norm(D1-D,'fro');
    e1(i) = norm(D1-D-h(i)*dD,'fro');
end
disp('Zeroth order Taylor error - should be ~1');
disp(median(log10(e0(2:end)./e0(1:end-1))));
disp('First order Taylor error - should be ~2');
disp(median(log10(e1(2:end)./e1(1:end-1))));

% Adjoint test
Z = randn(size(dD)) + 1i*randn(size(dD));
t = real( vec(Z)'*vec(dD) );
Jadj = PDEfunc_25D(PDEopts.JACOB_ADJ,vc,Q,Z,[],model,opts);
s = real(vec(Jadj)'*vec(dv));
disp('Adjoint test difference');
disp(abs(s-t)/abs(s));
