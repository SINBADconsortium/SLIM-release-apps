function [ A ] = funWindow1DavgBck( n, p, h )
%funWindow1DavgBck inverse avarage windowing sparse array for CARP-CG method
%
%   [ A ] = funWindow1DavgBck( N, P, H )
%
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      A = 'inverse' operator

    [ m os ys xs ] = pSPOT.pWindow.funWindow1DShape( n, p, h );

    A=sparse(n,m);
    for i=0:ys(1,2)-1
        A(ys(1,1)+i,ys(1,1)+i)=pSPOT.pWindow.Average1Dcr(i+1,ys(1,2),h);
    end
    for w=2:p-1
        for i=0:ys(w,2)-1
            A(ys(w,1)+i-2*h*(w-1),ys(w,1)+i)=pSPOT.pWindow.Average1Dmm(i+1,ys(w,2),h);
        end
    end
    for i=0:ys(p,2)-1
        A(ys(p,1)+i-2*h*(p-1),ys(p,1)+i)=pSPOT.pWindow.Average1Dcr(ys(p,2)-i,ys(p,2),h);
    end

    szs=size(A);
    assert(szs(1)==n,'funWindow1DavgBck: failed building the array');
    assert(szs(2)==m,'funWindow1DavgBck: failed building the array');

end
