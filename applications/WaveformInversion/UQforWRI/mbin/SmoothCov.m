function C = SmoothCov(n,a,b,c,dmax)

%% function to Generate the smooth CovMatrix:
%  C(i,j) = a * exp{- (r_i-r_j)^2 / 2b^2} + c * delta_{i,j}

if nargin < 5
  dmax = 3;
end

N = prod(n);

C = sparse(N,N);

for k = -dmax : dmax
  for l = 1 : dmax
    dxy  = norm([k,l])^2;
    mxy  = a * exp(-dxy/2/b^2);

    mvct = ones(N-k-abs(l)*n(1),1)* mxy;
    if k > 0
      idx = [n(1)-k+1:n(1)];
      for i = 1:n(2)-abs(l)-1
        mvct(idx) = 0;
        idx       = idx +n(1);
      end


    else
      if k < 0
        idx = [1:abs(k)];
        for i = 1:n(2)-abs(l)+1
          mvct(idx) = 0;
          idx       = idx +n(1);
        end
      end
    end

    nzero = k + abs(l)*n(1);
    mvct  = [zeros(nzero,1); mvct(:)];

    C     = C + spdiags(mvct,k+abs(l)*n(1),N,N);
  end
end

for k = 1 : dmax
    l = 0;
    dxy  = norm([k,l])^2;
    mxy  = a * exp(-dxy/2/b^2);

    mvct = ones(N-abs(k)-abs(l)*n(1),1)* mxy;
    if k > 0
      idx = [n(1)-k+1:n(1)];
      for i = 1:n(2)-abs(l)-1
        mvct(idx) = 0;
        idx       = idx +n(1);
      end


    else
      if k < 0
        idx = [1:abs(k)];
        for i = 1:n(2)-abs(l)+1
          mvct(idx) = 0;
          idx       = idx +n(1);
        end
      end
    end

    nzero = abs(k) + abs(l)*n(1);
    mvct  = [zeros(nzero,1); mvct(:)];

    C     = C + spdiags(mvct,k+abs(l)*n(1),N,N);

end


C = C'+C;

if length(c) == 1
  mvct = (a + c) * ones(N,1);
else
  mvct = a * ones(N,1) + c(:);
end
C    = C + spdiags(mvct,0,N,N);



end

function y = sub2ind_local(n,j,l)

    y = j + (l-1)*n(1);

end
