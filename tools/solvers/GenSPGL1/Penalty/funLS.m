function [f g] = funLS(r, params)
f = norm(r,2);
g = r./f;

end