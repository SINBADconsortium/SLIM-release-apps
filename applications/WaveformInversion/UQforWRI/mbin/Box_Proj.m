function m = Box_Proj(m,mmin,mmax)

idx = find(m<mmin);
m(idx) = mmin;
idx = find(m > mmax);
m(idx) = mmax;
