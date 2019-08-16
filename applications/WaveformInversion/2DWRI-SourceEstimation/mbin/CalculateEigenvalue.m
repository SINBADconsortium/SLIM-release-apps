function [lambda, A ,Pr] = CalculateEigenvalue(model,m)
    
    ot = model.o-model.nb.*model.d;
    dt = model.d;
    nt = model.n+2*model.nb;
    [zt,xt] = odn2grid(ot,dt,nt);
    Nt      = prod(nt);
    freq    = model.freq;

    % data size
    nrec   = length(model.zrec)*length(model.xrec);
    nfreq  = length(model.freq);
    Pr     = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
    Pr     = sparsedouble(Pr);
    Px     = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
    
    mx = Px*m;
    
    [fm,df,ddf]    = input2helm_param(mx,freq);
    [Hk,A,B]       = Helm2D(fm,ot,dt,nt,model.nb);
    
    L      = Hk'\Pr';
    lambda = svds(L,1);
    
        
    
    
