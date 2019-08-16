function [ A ] = funWindow1DfdFor( n, p, h )
%funWindow1DfdFor forward swap windowing sparse array for finite difference algorithms
%
%   [ A ] = funWindow1DfdFor( N, P, H )
%
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      A = 'inverse' operator

    [ m os ys xs ] = pSPOT.pWindow.funWindow1DShape( n, p, h );

    A=sparse(m,n);
    for w=1:p
        for i=0:ys(w,2)-1
            A(ys(w,1)+i,ys(w,1)+i-2*h*(w-1))=1;
        end
    end

    szs=size(A);
    assert(szs(1)==m,'funWindow1DfdFor: failed building the array');
    assert(szs(2)==n,'funWindow1DfdFor: failed building the array');

end
