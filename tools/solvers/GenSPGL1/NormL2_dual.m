function d = NormL2_dual(x,weights)

d = norm(x./weights,2);
