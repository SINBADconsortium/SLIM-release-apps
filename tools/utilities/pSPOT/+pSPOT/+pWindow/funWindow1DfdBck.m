function [ A ] = funWindow1DfdBck( n, p, h )
%funWindow1DfdBck inverse swap windowing sparse array for finite difference algorithms
%
%   [ A ] = funWindow1DfdBck( N, P, H )
%
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      A = 'inverse' operator

    [ m os ys xs ] = pSPOT.pWindow.funWindow1DShape( n, p, h );

    A=sparse(n,m);
    for w=1:p
        for i=0:xs(w,2)-1
            A(xs(w,1)+i,xs(w,1)+i+2*h*(w-1))=1;
        end
    end

    szs=size(A);
    assert(szs(1)==n,'funWindow1DfdBck: failed building the array');
    assert(szs(2)==m,'funWindow1DfdBck: failed building the array');

end
