function c = prx(m)

if nargin < 1
    m = size(get(gcf,'colormap'),1);
end

if (mod(m,2) == 0)
    m1 = m*0.5;
    taper = (0:m1-1)'/max(m1-1,1);
    taper2 = (2*taper - 1);
    zeropart = (taper2 < 0);
    taper2(zeropart) = zeros(sum(zeropart),1);
    
    r = [taper; ones(m1,1)];
    g = [taper; flipud(taper)];
    b = [taper; flipud(taper2)];
else
    m1 = floor(m*0.5);
    taper = (0:m1-1)'/max(m1,1);
    taper2 = (2*taper - 1);
    zeropart = (taper2 < 0);
    taper2(zeropart) = zeros(sum(zeropart),1);

    r = [taper; ones(m1+1,1)];
    g = [taper; 1; flipud(taper)];
    b = [taper; 1; flipud(taper2)];
end

c = [r g b];

end

