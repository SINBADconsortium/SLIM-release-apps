function vs = vsmooth(v,d)

n = size(v);

p1 = spdiags(repmat([.25 .5 .25],n(1),1),-1:1,n(1),n(1)); 
p1(1,1) = p1(1,1)+.25; p1(n(1),n(1)) = p1(n(1),n(1))+.25;
p2 = spdiags(repmat([.25 .5 .25],n(2),1),-1:1,n(2),n(2)); 
p2(1,1) = p2(1,1)+.25; p2(n(2),n(2)) = p2(n(2),n(2))+.25;

vs = (p1^d)*v*(p2^d);