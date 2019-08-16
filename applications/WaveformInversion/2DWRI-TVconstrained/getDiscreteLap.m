function L = getDiscreteLap(n,h)

D1 = spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1)); D1(1,1) = -1; D1(end,end) = -1;
D2 = spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2)); D2(1,1) = -1; D2(end,end) = -1;
L = (kron(speye(n(2)),D1) + kron(D2,speye(n(1))))/h/h;