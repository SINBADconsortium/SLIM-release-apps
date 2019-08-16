function [model,pm,Ps,Pr,q,d,ssW] = setup_problem(example)
%load or define data
%also define example specific algorithm parameters in pm and model info in model

switch example
    
    case {'blob'}
        %% blob example
        model.example = example;
        %model parameters and data
        h = 10; model.h = h; %mesh size in meters
        n = round([203 201]); model.n = n; %(z,x)
        N = prod(n); model.N = N;
        model.zt = 0:h:h*(n(1)-1);
        model.xt = 0:h:h*(n(2)-1);
        %define true velocity model
        vbmin = 1500;
        vbmax = 3000;
        vgrad = (vbmin:(vbmax-vbmin)/(n(1)-1):vbmax)'*ones(1,n(2));
        [xx,zz] = meshgrid(1:model.n(2),1:model.n(1));
        dv = 0*xx;
        cx = round(n(2)/2);
        cz = round(5*n(1)/8);
        r = 300/h;
        dv((xx-cx).^2 + (zz-cz).^2 <= r^2) = 4500;
        vtrue = vgrad;
        vtrue(dv>0) = dv(dv>0); %add high velocity blob to vgrad
        model.vtrue = vtrue;
        %define initial velocity to be background
        model.vinit = vgrad;
        model.vmin = 1000*ones(n(1),n(2));
        model.vmax = 5000*ones(n(1),n(2));
        model.mtrue = 1./model.vtrue.^2;
        model.minit = 1./model.vinit.^2;
        model.mmin = 1./model.vmax.^2;
        model.mmax = 1./model.vmin.^2;
        %define frequencies (Hz) to be used
        model.freq = 3:1:33;
        nf = length(model.freq); model.nf = nf;
        model.Vb = [(1:nf-1)' (2:nf)']; %every row defines a frequency batch
        model.batches = size(model.Vb,1);
        model.f0 = 30; %peak freq of ricker wavelet
        model.t0 = 0; %phase shift of wavelet in seconds
        %receivers and sources near the top
        model.zsrc = 5*h:5*h:195*h;
        ns = length(model.zsrc); model.ns = ns;
        model.xsrc = 20*h*ones(1,ns);
        model.zrec = 5*h:2*h:195*h;
        nr = length(model.zrec); model.nr = nr;
        model.xrec = 180*h*ones(1,nr);
        %Ricker wavelet (freq domain) with peak frequency f0, phase shift t0
        q = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*exp(1i*2*pi*model.freq*model.t0);
        model.vweights = ones(1,nf); %frequency weights (not important for small freq batches)
        %define weights for simultaneous shots
        model.nsim = ns; %number of simultaneous shots
        model.redraw = 0; %redraw flag: nsim random simultaneous shots if 1, ssW stays same if 0
        model.rngshots = rng('shuffle'); %can also load some seed here
        rng(model.rngshots);
        %ssW = randn(ns,model.nsim); %shot weights
        ssW = speye(ns,model.nsim); %sequential shots (also set nsim = ns and redraw=0 in this case)
        % sampling operators
        Ps = sparse((1:ns)',round(model.zsrc/model.h+1) + ...
            round(model.xsrc/model.h)*n(1),ones(ns,1),ns,N);
        Pr = sparse((1:nr)',round(model.zrec/model.h+1) + ...
            round(model.xrec/model.h)*n(1),ones(nr,1),nr,N);
        %boundary condition parameters
        Xint = [0 ones(1,n(1)-2) 0];
        Xint = Xint(:)*[0 ones(1,n(2)-2) 0];
        model.Xint = Xint(:);
        Xbnd = (1-Xint);
        Xbnd(1,1) = Xbnd(1,1)*2; Xbnd(end,end) = Xbnd(end,end)*2;
        Xbnd(end,1) = Xbnd(end,1)*2; Xbnd(1,end) = Xbnd(1,end)*2;
        model.Xbnd = Xbnd(:)/h;
        %generate data
        L = getDiscreteLap(n,h); model.L = L;
        model.rngnoise = rng('shuffle'); %can load some seed here
        model.dsig = .05; %try .05 for a noisy data test
        %Define mass matrix as function of frequency and model including Robin BC
        model.M = @(ff,mm) (2*pi*ff)^2*spdiags(model.Xint(:).*mm(:),0,N,N) - ...
            1i*(2*pi*ff)*spdiags(model.Xbnd(:).*sqrt(mm(:)),0,N,N);
        disp('generating synthetic data...');
        d = generateData(model,q,Ps,Pr); %generate data
        
        %algorithm parameters
        %pm.dpw = (1./((1500:(2500/(n(1)-1)):4000)'/(1500))).^2; %decreasing depth weights (max must be 1)
        pm.dpw = ones(n(1),1); %uniform depth weights
        [D,E] = getDiscreteGrad(n(1),n(2),h,h,pm.dpw); %weights are built in to D
        model.D = D; model.E = E; %discrete grad and forward differences
        pm.TV = @(t) sum(sqrt(E'*((D*t(:)).^2))); %function to compute TV
        pm.mu = 1; %numerical parameter (doesn't change objective)
        pm.lambda = 1; %set penalty parameter for PDE misfit
        pm.tau = pm.TV(model.mtrue); %total variation of true model
        pm.tau = pm.tau*.875; %set parameter for TV constraint
        %pm.tau = pm.tau*1000; %large enough so that TV penalty has no effect
        %outer trust region/step size parameter c and how to adjust it dynamically
        pm.cmin = 1e-14; pm.cmax = 1e16;
        pm.c1 = 2; %factor to decrease c if objective decreasing enough
        pm.c2 = 10; %factor to increase c if objective doesn't decrease enough
        pm.sigma = .1; %what fraction of ideal objective decrease is acceptable
        %convex subproblem parameters
        pm.itol = 1e-4; %tolerance for inner iteration stopping condition
        pm.miniits = 10; %minimum number of inner iterations
        pm.maxiits = 2000; %maximum number of allowed inner iterations
        pm.admax = h^2/8; %conservative lower bound for 1/||D'*D||
        %SGP parameters
        pm.otol = 1e-4; %tolerance for outer iteration stopping condition
        pm.maxoits = 25; %maximum number of allowed outer iterations
        pm.minoits = min(2,pm.maxoits); %minimum number of outer iterations
        
    
    case {'BPC'} 
        %% BP center example
        model.example = example;
        %model parameters and data 
        load BPCtrue 
        h = 20; model.h = h; %mesh size in meters
        n = size(BPCtrue); model.n = n; %(z,x)
        N = prod(n); model.N = N;
        model.zt = 0:h:h*(n(1)-1);
        model.xt = 0:h:h*(n(2)-1);
        %set up velocity model
        model.vtrue = BPCtrue;
        %define initial velocity by smoothing model.vtrue
        model.vinit = vsmooth(model.vtrue,172);
        model.vmin = 1400*ones(n(1),n(2));
        model.vmax = 5000*ones(n(1),n(2));
        model.mtrue = 1./model.vtrue.^2;
        model.minit = 1./model.vinit.^2;
        model.mmin = 1./model.vmax.^2;
        model.mmax = 1./model.vmin.^2;
        %define frequencies (Hz) to be used
        model.freq = 3:1:20;
        nf = length(model.freq); model.nf = nf;
        model.Vb = [(1:nf-1)' (2:nf)']; %every row defines a frequency batch
        model.batches = size(model.Vb,1);
        model.f0 = 15; %peak freq of ricker wavelet
        model.t0 = 0; %phase shift of wavelet in seconds
        %receivers and sources near the top
        model.xsrc = 50*h:4*h:(n(2)-50)*h;
        ns = length(model.xsrc); model.ns = ns;
        model.zsrc = 2*h*ones(1,ns);
        model.xrec = 2*h:2*h:(n(2)-2)*h;
        nr = length(model.xrec); model.nr = nr;
        model.zrec = 3*h*ones(1,nr);
        %Ricker wavelet (freq domain) with peak frequency f0, phase shift t0
        q = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*exp(1i*2*pi*model.freq*model.t0);
        model.vweights = ones(1,nf); %frequency weights (not important for small freq batches)
        %define weights for simultaneous shots
        model.nsim = 2; %number of simultaneous shots
        model.redraw = 1; %redraw flag: nsim random simultaneous shots if 1, ssW stays same if 0
        model.rngshots = rng('shuffle'); %can also load some seed here
        rng(model.rngshots);
        ssW = randn(ns,model.nsim); %shot weights
        %ssW = speye(ns,model.nsim); %sequential shots (also set nsim = ns and redraw=0 in this case)
        % sampling operators
        Ps = sparse((1:ns)',round(model.zsrc/model.h+1) + ...
            round(model.xsrc/model.h)*n(1),ones(ns,1),ns,N);
        Pr = sparse((1:nr)',round(model.zrec/model.h+1) + ...
            round(model.xrec/model.h)*n(1),ones(nr,1),nr,N);
        %boundary condition parameters
        Xint = [0 ones(1,n(1)-2) 0];
        Xint = Xint(:)*[0 ones(1,n(2)-2) 0];
        model.Xint = Xint(:);
        Xbnd = (1-Xint);
        Xbnd(1,1) = Xbnd(1,1)*2; Xbnd(end,end) = Xbnd(end,end)*2;
        Xbnd(end,1) = Xbnd(end,1)*2; Xbnd(1,end) = Xbnd(1,end)*2;
        model.Xbnd = Xbnd(:)/h;
        %generate data
        L = getDiscreteLap(n,h); model.L = L;
        model.rngnoise = rng('shuffle'); %can load some seed here
        model.dsig = 0; %try .05 for a noisy data test
        %Define mass matrix as function of frequency and model including Robin BC
        model.M = @(ff,mm) (2*pi*ff)^2*spdiags(model.Xint(:).*mm(:),0,N,N) - ...
            1i*(2*pi*ff)*spdiags(model.Xbnd(:).*sqrt(mm(:)),0,N,N);
        disp('generating synthetic data...');
        d = generateData(model,q,Ps,Pr); %generate data
        
        %algorithm parameters
        %pm.dpw = (1./((1500:(2500/(n(1)-1)):4000)'/(1500))).^2; %decreasing depth weights (max must be 1)
        pm.dpw = ones(n(1),1); %uniform depth weights
        [D,E] = getDiscreteGrad(n(1),n(2),h,pm.dpw); %weights are built in to D
        model.D = D; model.E = E; %discrete grad and forward differences
        pm.TV = @(t) sum(sqrt(E'*((D*t(:)).^2))); %function to compute TV
        pm.mu = 1; %numerical parameter (doesn't change objective)
        pm.lambda = 1; %set penalty parameter for PDE misfit
        pm.tau = pm.TV(model.mtrue); %total variation of true model
        pm.tau = pm.tau*.9; %set parameter for TV constraint
        %pm.tau = pm.tau*1000; %large enough so that TV penalty has no effect
        %outer trust region/step size parameter c and how to adjust it dynamically
        pm.cmin = 1e-14; pm.cmax = 1e16;
        pm.c1 = 2; %factor to decrease c if objective decreasing enough
        pm.c2 = 10; %factor to increase c if objective doesn't decrease enough
        pm.sigma = .1; %what fraction of ideal objective decrease is acceptable
        %convex subproblem parameters
        pm.itol = 1e-4; %tolerance for inner iteration stopping condition
        pm.miniits = 10; %minimum number of inner iterations
        pm.maxiits = 2000; %maximum number of allowed inner iterations
        pm.admax = h^2/8; %conservative lower bound for 1/||D'*D||
        %SGP parameters
        pm.otol = 1e-4; %tolerance for outer iteration stopping condition
        pm.maxoits = 25; %maximum number of allowed outer iterations
        pm.minoits = min(2,pm.maxoits); %minimum number of outer iterations
         
        
    otherwise
        disp([example ' not defined in setup_problem.m']);
        model = []; pm = []; Ps = []; Pr = []; q = []; d = []; ssW = [];
end