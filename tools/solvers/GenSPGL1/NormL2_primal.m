function p = NormL2_primal(x,weights)

p = norm(x.*weights,2);
