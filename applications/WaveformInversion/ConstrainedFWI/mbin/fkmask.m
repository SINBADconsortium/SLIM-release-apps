function maskvector = fkmask (M, N, lbound, ubound, aspect)
    dims = [M N];
    MEdims = 2*dims;
    
    bx = (MEdims(1)) / 2;
    bz = (MEdims(2)) / 2;
    
    [x, z] = meshgrid(-bz:bz-1, -bx:bx-1);
    x = aspect * (x / bx);
    z = z / bz;
    
    dist = sqrt(x .^2 + z .^2);
    
    mask = cos((pi / 2).*(dist - lbound)./(ubound - lbound)).^2;
    
    mask(dist >= ubound) = 0;
    mask(dist <= lbound) = 1;
    
    maskvector = vec(mask);
end
