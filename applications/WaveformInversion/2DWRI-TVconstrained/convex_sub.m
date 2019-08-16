function u = convex_sub(model,pm,m,g,Hc)
%PDHG variant for solving
%min_u <u,g> + .5<u,Hc u> s.t. m_ij + u_ij in [b_ij,B_ij] and ||m+u||TV <= tau 

%INPUT:
%m = current model estimate
%g = gradient
%H = Hessian approximation
%c = parameter for damping H (Hc = H + 2cI)

%OUTPUT:
%estimated model increment u

%initialize stuff
%u = -g./Hc; %unconstrained Gauss Newton update
%return
D = model.D;
E = model.E;
u = -g./Hc; 
p = zeros(size(D,1),1);
fac = max(Hc);

%parameters
a = 1/fac;
d = pm.admax*fac;
tau = pm.tau;

%track objective and constraint
%obj = g'*u + .5*u'*(Hc.*u); 
%con = model.TV(m(:)+u(:)) - tau;

%iterate
itol = pm.itol;
ires = 2*itol;
urel = 2*itol;
prel = 2*itol;
iit = 1;
maxiit = pm.maxiits;
miniit = min(maxiit,pm.miniits);
while ( iit <= maxiit && (ires>itol || iit<=miniit) )
    
    %p update
    t = p + d*D*(m(:) + u);
    p_prev = p;
    p = t - proj_E(E,t,tau*d); %vector l1 ball projection
    
    %u update (includes bound constraints)
    u_prev = u;
    u = max(model.mmin(:)-m(:),min(model.mmax(:)-m(:),...
        (-g + u/a - D'*(2*p - p_prev))./(Hc+1/a)));

    %stop when u not changing much and when constraint approx satisfied
    %obj(iit+1) = g'*u + .5*u'*(Hc.*u); 
    urel(iit+1) = norm(u-u_prev)/(norm(u)+eps);
    prel(iit+1) = norm(p-p_prev)/(norm(p)+eps); %check how much p is changing
    %con(iit+1) = model.TV(m(:)+u(:)) - tau;
    ires = max(urel(iit+1),prel(iit+1));
    
    %option to continue iterating (this is mainly for debugging, can remove later)
    if (iit==maxiit+1) %change to (iit==maxiits) to plot intermediate progress
        %plot stuff
        %figure(90)
        %clf;
        %plot(obj); title('<g,u> + .5*<u,Hc*u>');
        figure(91)
        clf;
        plot(urel); title('relative change in primal')      
        figure(92)
        clf;
        plot(prel); title('|relative change in dual');
        
        %option to keep iterating when maxiit reached
        mt = input('hit enter to finish or enter number of inner its to continue: ');
        if (~isempty(mt))
            maxiit = maxiit + mt;
        end
    end
    iit = iit + 1;
end
if (iit-1 == maxiit)
    disp(['stopping condition not reached in ' num2str(iit-1) ' inner iterations, ires = ' num2str(ires)]);
else
    disp(['stopping condition was reached in ' num2str(iit-1) ' inner iterations, ires = ' num2str(ires)]);
end
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
