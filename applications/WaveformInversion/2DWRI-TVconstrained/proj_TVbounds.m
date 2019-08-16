function m = proj_TVbounds(m0,tau,h,dpw,b,B)
%project m0 onto TV ball of radius tau
%intersected with bound constraints m_ij in [b_ij,B_ij],
%where h is the mesh width 
%and dpw are depth weights incorporated into the definition of TV
%by solving min_m .5||m - m0||^2 s.t. ||m||_TV <= tau, m_ij in [b_ij,B_ij]
%Note this is very similar to convex_sub

%get discrete gradient
[n1,n2] = size(m0);
[D,E] = getDiscreteGrad(n1,n2,h,h,dpw); %D includes depth weights dpw

%algorithm parameters
fac = 1/h/h;
a = 1*fac;
d = h^2/8/fac; %must have a*d < 1/||D'*D||

%initialize stuff
N = n1*n2;
u = zeros(N,1);
p = zeros(size(D,1),1);

%iterate
itol = 1e-5;
ires = 2*itol;
urel = 2*itol;
prel = 2*itol;
iit = 1;
maxiit = 5000;
miniit = 10;
while ( iit <= maxiit && (ires>itol || iit<=miniit) )
    
    %p update
    t = p + d*D*u;
    p_prev = p;
    p = t - proj_E(E,t,tau*d); %vector l1 ball projection
    
    %u update
    u_prev = u;
    u = max(b(:),min(B(:),(a*m0(:) + u - a*D'*(2*p - p_prev))./(a+1)));
    
    %stop when u and p not changing much
    urel(iit+1) = norm(u-u_prev)/(norm(u)+eps);
    prel(iit+1) = norm(p-p_prev)/(norm(p)+eps);
    ires = max(urel(iit+1),prel(iit+1));

    iit = iit + 1;
end
disp([num2str(iit-1) ' totel iterations, ires = ' num2str(ires)]);
m = reshape(u,n1,n2);

end


function px = proj_E(E,x,ep)
%compute orthogonal projection onto vector l1 ball
%defined in terms of E by sum(sqrt(E'*(x.^2))) <= ep

%get magnitudes
xm = sqrt(E'*(x.^2));

%done if inside region, otherwise project magnitudes onto simplex
if (sum(xm) <= ep)
    px = x;
else
    %normalize the vectors that are grouped by E
    xn = x./(E*xm+eps);
    
    %project xm onto simplex defined by sum to c constraint
    pxm = proj_simplex(xm,ep);
    
    %combine projection with original directions to get px
    px = (E*pxm).*xn;
end
end


function px = proj_simplex(x,ep)
%using bisection to figure out how far perpendicular to hyperplane to move
%so that nonnegative projection sums to ep
Tmin = min(x);
Tmax = max(x);
Tmid = (Tmin+Tmax)/2;
px = max(0,x-Tmid);
err = sum(px) - ep;
tol = 1e-12;
maxit = 50; %should be plenty
it = 1;
while (it <= maxit && abs(err) > tol)
    if (err > 0)
        Tmin = Tmid;
        Tmid = (Tmid+Tmax)/2;
    else
        Tmax = Tmid;
        Tmid = (Tmin+Tmid)/2;
    end
    px = max(0,x-Tmid);
    err = sum(px) - ep; %could probably do this more efficiently
    it = it + 1;
end
end
