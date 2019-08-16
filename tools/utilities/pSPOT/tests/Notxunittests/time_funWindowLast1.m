function test_funWindowLast1()
    hh=[1 1 3];
    mm=[100 100 1000];
    tt=[1 2];
    p=parpool_size();

    fprintf('%3s | %4s %4s %4s %4s | %s |','p','h','m','l','N','t');
    fprintf('%6s %6s %5s %6s %5s %6s %5s\n','td','tc','tc/td','to','to/td','tk','tk/td')
    for h=hh(1):hh(2):hh(3)
        for m=mm(1):mm(2):mm(3)
            for t=tt(1):tt(2)
                % prepare input vector
                l=m+2*h*p;
                odims=[m,l];

                % make x
                x=distributed.randn(odims);
                x=x(:);
                n=prod(size(x));
        
                % warmup smpd
                spmd
                    xx=pSPOT.pWindow.funWindowLast1HaloMakeCodist(x(:),l,p,h);
                    xx=pSPOT.pWindow.funWindowLast1HaloAverageCodist(xx,l,p,h);
                    x=pSPOT.pWindow.funWindowLast1HaloDropCodist(xx,l,p,h);
                end

                % run function form with MPI for codistributed
                [ N xs ys ] = pSPOT.pWindow.funWindowLast1HaloShape(l,p,h);
                spmd
                    tic;
                    yl=pSPOT.pWindow.funWindowLast1HaloMakeCodist(x,l,p,h);
                    zl=pSPOT.pWindow.funWindowLast1HaloAverageCodist(yl,l,p,h);
                    vl=pSPOT.pWindow.funWindowLast1HaloDropCodist(zl,l,p,h);
                    tc=toc;
                end
                tc=max(cell2mat(SeisDataContainer.utils.Composite2Cell(tc)));

                % run function form with MPI for distributed
                [ N xs ys ] = pSPOT.pWindow.funWindowLast1HaloShape(l,p,h);
                tic;
                y=pSPOT.pWindow.funWindowLast1HaloMakeDist(x,l,p,h);
                z=pSPOT.pWindow.funWindowLast1HaloAverageDist(y,l,p,h);
                v=pSPOT.pWindow.funWindowLast1HaloDropDist(z,l,p,h);
                td=toc;
        
                % run operator form with MPI for distributed
                H=opdWindowLast1Halo(n,l,p,h);
                E=opdWindowLast1HaloAverage(n,l,p,h);
                tic;
                Y=H*x;
                Z=E*Y;
                V=H'*Z;
                to=toc;
        
                % run operator form with Kronecker for distributed
                A=oplWindow1Davg(l,p,h);
                D=opDirac(m);
                W=oppKron2Lo(A,D);
                tic;
                YY=W*x;
                ZZ=W*W'*YY;
                VV=W'*ZZ;
                tk=toc;
        
                fprintf('%3d | %4d %4d %4d %4d | %d |',p,h,m,l,N,t);
                fprintf('%6.3f %6.3f %5.2f %6.3f %5.2f %6.3f %5.2f\n',td,tc,tc/td,to,to/td,tk,tk/td)
            end
        end
    end
end
