function p = NewtonMcMC(fh0, fq0, fh1, fq1)

    fh0 = -fh0;
    fh1 = -fh1;
    fq0 = -fq0;
    fq1 = -fq1;
    
    dfh = fh1 - fh0
    dfq = fq0 - fq1
    
    p   = min(1, exp(dfh+dfq));
