function [f g] = hubers(x, thresh, slope)


if(nargin < 2)
    thresh = 1;
    slope = 1;
end

if(nargin < 3)
   slope = 1;  
end

x = x(:);

Ia = find(abs(x) <= thresh);
Ib = find(abs(x) >  thresh);

f = 0;
if ~isempty(Ia)
    f = .5*x(Ia)'*x(Ia)*slope/thresh;
end
if ~isempty(Ib)
    f = f + slope*sum(abs(x(Ib)) - .5*thresh);
end
f = gather(f);

g     = x;
g(Ia) = slope*x(Ia)/thresh;
g(Ib) = slope*sign(x(Ib));


end
