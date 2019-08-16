function test_suite = test_ToolsAlgorithms3DFreqModeling
% Unit tests for PDEfunc, 3d mode
% 
%
% Curt Da Silva, 2015 
%
initTestSuite;

function model = setup
    nx = 50; ny = 50; nz = 50;
    nsrc = 1; nrec = 100;
    model.o = [0 0 0];
    model.n = [nx ny nz];
    model.d = [10 10 10];
    dx = model.d(1); dy = model.d(2); dz = model.d(3);
    [x,y,z] = odn2grid(model.o,model.d,model.n);
    
    model.freq = 4;
    model.f0 = 10;
    model.t0 = 0;
    model.unit = 'm/s';
    model.zsrc = 100;
    model.xsrc = linspace(4*dx,max(x)-4*dx,nsrc);
    model.ysrc = linspace(4*dy,max(y)-4*dy,nsrc);
    model.zrec = 100;
    model.xrec = linspace(4*dx,max(x)-4*dx,nrec);
    model.yrec = linspace(4*dy,max(y)-4*dy,nrec);
    model.v = 2000*ones(model.n);
    lsopts = LinSolveOpts();
    lsopts.tol = 1e-10;
    lsopts.maxit = 2000;
    lsopts.maxinnerit = 5;
    lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
    lsopts.precond = LinSolveOpts.PREC_MLGMRES;
    
    opts = struct;
    opts.free_surface = false;
    opts.disp_output = false;
    opts.scheme = PDEopts.HELM3D_OPERTO27;
    opts.solve_opts = lsopts;
    opts.pml = 10;
    pdeopts = PDEopts;
    pdeopts.helm_pml = 10;
    opts.pdefunopts = pdeopts;
    opts.lsopts = lsopts;
    model.opts = opts;
    
   
function testdH(model)
   tol = 0.25;
   v = model.v; 
   
   if strcmp(model.unit,'s2/m2')
       v = (v.^(-2));
   end
   freq = model.freq;
      
   dv = (max(vec(v))/10)*randn(size(v));   
      
   opts = model.opts;
   
   [H,comp_grid,dH,ddH] = discrete_helmholtz(v,model,freq(1),opts);
   dH = dH*opDiag(comp_grid.phys2comp*vec(dv));
   ddH = ddH*opDiag((comp_grid.phys2comp*vec(dv)).^2);
   h = 10.^(-6:-1);
   e0 = zeros(length(h),1); e1 = zeros(length(h),1); e2 = zeros(length(h),1);
   test_vec = randn(size(H,2),1);
   for i=1:length(h)
      v1 = v+h(i)*dv;           
      H1 = discrete_helmholtz(v1,model,freq,opts);
      e0(i) = norm(H1*test_vec-H*test_vec);
      e1(i) = norm(H1*test_vec-H*test_vec-h(i)* dH*test_vec);      
      e2(i) = norm(H1*test_vec-H*test_vec-h(i)* dH*test_vec -h(i)^2/2*ddH*test_vec);      
   end
   h0 = median(log10(e0(2:end)./e0(1:end-1)));   
   assert(abs(h0-1) < tol);
   switch model.unit
     case 'm/s'
          h1 = median(log10(e1(2:end)./e1(1:end-1)));
          assert(abs(h1-2)<tol);
          h2 = max(log10(e2(2:end)./e2(1:end-1)));
          assert(abs(h2-3)<tol);
     case 's2/m2'
       % helmholtz is linear in slowness squared
          assert(max(e1)<1e-8);
   end
   
   
function testdU(model)
   tol = 0.3;
   v = model.v; 
   if strcmp(model.unit,'s2/m2')
       v = (v.^(-2));
   end      
   
   freq = model.freq;
      
   dv = (max(vec(v))/10)*randn(size(v));   
      
   opts = model.opts;
      
   [H,comp_grid,dH] = discrete_helmholtz(v,model,freq(1),opts);
   dHdv = dH*opDiag(comp_grid.phys2comp*vec(dv));
   Q = zeros(size(v)); Q(round(model.n(1)/2),round(model.n(2)/2),round(model.n(3)/2)) = 1;
   Q = comp_grid.phys2comp*vec(Q);
   U = H\Q;
   dU = H\( dHdv*(-U) );
   h = 10.^(-3:0);
   e0 = zeros(length(h),1); e1 = zeros(length(h),1);
   for i=1:length(h)
       v1 = v+h(i)*dv;
       assert(min(vec(v1)) > 0);
       H1 = discrete_helmholtz(v1,model,freq(1),opts);
       U1 = H1\Q;
       e0(i) = norm(U1-U);
       e1(i) = norm(U1-U-h(i)*dU);
   end
   h0 = median(log10(e0(2:end)./e0(1:end-1)));
   h1 = median(log10(e1(2:end)./e1(1:end-1)));
   
   assert(abs(h0-1) < tol);
   assert(abs(h1-2) < tol);

   
function testJacobian(model)
   tol = 0.3;
   v = model.v; freq = model.freq; 
   if strcmp(model.unit,'s2/m2')
       v = (v.^(-2));
   end      
   
   dv = (max(vec(v))/100)*randn(size(v)); 
   dv([1 end],:,:) = 0; dv(:,[1 end],:) = 0; dv(:,:,[1 end]) = 0;
   opts  = model.opts;
      
   nsx = length(model.xsrc); nsy = length(model.ysrc);
   nsrc = nsx*nsy;
   nrx = length(model.xrec); nry = length(model.yrec);
   nrec = nrx*nry;
   nfreq = length(model.freq);
   Q = speye(nsrc);
   
   D0 = PDEfunc(PDEopts.FORW_MODEL,v,Q,[],[],model,opts);
   dD = PDEfunc(PDEopts.JACOB_FORW,v,Q,dv,[],model,opts);
   h = 10.^(-3:0);
   e0 = zeros(length(h),1); e1 = zeros(length(h),1);
   for i=1:length(h)
       v1 = v+h(i)*dv;
       assert(min(vec(v1)) > 0);
       D1 = PDEfunc(PDEopts.FORW_MODEL,v1,Q,[],[],model,opts);
       e0(i) = norm(vec(D1-D0));
       e1(i) = norm(vec(D1-D0-h(i)*dD));
   end
   h0 = median(log10(e0(2:end)./e0(1:end-1)));
   h1 = median(log10(e1(2:end)./e1(1:end-1)));
   
   assert(abs(h0-1) < tol);
   assert(abs(h1-2) < tol);

   % Test adjoint
   Z = randn(size(dD))+1i*randn(size(dD));
   s = real( vec(dD)'*vec(Z) );
   y = PDEfunc(PDEopts.JACOB_ADJ,v,Q,Z,[],model,opts);
   t = vec(dv)'*vec(y);
   
   % Not completely accurate, since helmholtz matrices inverted inexactly
   assert(abs(s-t)<1e-2*max(abs(s),abs(t)));
   
function testObjective(model)
   tol = 0.25;
   v = model.v;   
   opts  = model.opts;  
   opts.solve_opts.tol = 1e-10;
   nsrc = length(model.xsrc)*length(model.ysrc);
   Q = speye(nsrc);
   
   dv = max(vec(v))/100*randn(model.n);
   dv([1 end],:,:) = 0; dv(:,[1 end],:) = 0; dv(:,:,[1 end]) = 0;
   Dobs = F3d(v+dv,Q,model,opts);
   dv = max(vec(v))/10*randn(model.n);
   dv([1 end],:,:) = 0; dv(:,[1 end],:) = 0; dv(:,:,[1 end]) = 0;
   if parpool_size()==0
       obj = @(h) PDEfunc(PDEopts.OBJ,v+h*dv,Q,[],Dobs,model,opts);
   else
       obj = @(h) PDEfunc_dist(PDEopts.OBJ,v+h*dv,Q,[],Dobs,model,opts);
   end
   [f0,g0] = obj(0);
         
   h = 10.^(-5:0);
   e0 = zeros(length(h),1); e1 = zeros(length(h),1); e2 = zeros(length(h),1);
   if parpool_size()==0
       Hdv = PDEfunc(PDEopts.HESS,v,Q,dv,Dobs,model,opts);
   else
       Hdv = PDEfunc_dist(PDEopts.HESS,v,Q,dv,Dobs,model,opts);
   end
   df = vec(g0)'*vec(dv);
   d2f = vec(dv)'*vec(Hdv);
   for i=1:length(h)
       v2 = v+h(i)*dv;
       assert(min(vec(v2)) > 0);
       f1 = obj(h(i));
       e0(i) = norm(f1-f0);
       e1(i) = norm(f1-f0-h(i)*df);
       e2(i) = norm(f1-f0-h(i)*df-h(i)^2/2*d2f);
   end
   h0 = median(diff(log10(e0)));
   h1 = median(diff(log10(e1)));
   h2 = median(diff(log10(e2)));
   assert(abs(h0-1) < tol || h0 > 1);
   assert(abs(h1-2) < tol || h1 > 2);
   assert(abs(h2-3) < tol || h2 > 3);         
   
  
   
function testGaussNewton(model)
   tol = 0.25;
   v = model.v; freq = model.freq;  
   if strcmp(model.unit,'s2/m2')
       v = (v.^(-2));
   end      
   
   dv = (max(vec(v))/100)*randn(size(v)); 
   dv([1 end],:,:) = 0; dv(:,[1 end],:) = 0; dv(:,:,[1 end]) = 0;
   opts  = model.opts;
   
   nsx = length(model.xsrc); nsy = length(model.ysrc);
   nsrc = nsx*nsy;
   nrx = length(model.xrec); nry = length(model.yrec);
   nrec = nrx*nry;
   nfreq = length(model.freq);
   Q = speye(nsrc);   
   dD = PDEfunc(PDEopts.JACOB_FORW,v,Q,dv,[],model,opts);
   Hdm = PDEfunc(PDEopts.JACOB_ADJ,v,Q,dD,[],model,opts);
   Hdm2 = PDEfunc(PDEopts.HESS_GN,v,Q,dv,[],model,opts);
   assert(norm(Hdm - Hdm2) < 1e-10*norm(Hdm));
   
