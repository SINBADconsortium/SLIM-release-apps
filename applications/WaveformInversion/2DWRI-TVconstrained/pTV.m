function [m,energy,oits,c] = pTV(model,pm,b,minit,Ps,Pr,q,d,ssW,cinit)
%min_{m,u} .5||Pu-d||^2 + .5lam^2||A(m)u-q||^2
%s.t. ||m||_TV <= tau and m_ij in [b_ij,B_ij]

%grab just the frequency data for current batch
freq_ind = model.Vb(b,:);
freq = model.freq(freq_ind);
q = q(freq_ind);
d = d(freq_ind);

%initialize stuff
Q = speye(model.ns); %source weight (simultaneous sources)
n = model.n;
m = minit;
c = max(pm.cmin,cinit); %if H=0, 1/2c would be gradient descent time step
if (model.redraw == 1)
    ssW = drawWeights(model);
end

%evaluate objective obj, define gradient g and Gauss Newton Hessian H
[obj,g,H] = computeGradient(m,model,pm,Q,ssW,Ps,Pr,freq_ind,freq,q,d);
Hc = H+2*c;
energy = obj; %track objective for plotting

%outer loop of Gauss-Newton-like updates
ores = 2*pm.otol;
oit = 1;
discarded_oits = 0;
maxoits = pm.maxoits;
minoit = min(pm.minoits,maxoits); %do at least 2 iterations
while ( oit <= maxoits && (ores>pm.otol || oit <= minoit) )
    
    %solve convex subproblem using modified PDHG, a linearized variant of ADMM
    %min_u <u,g> + .5<u,(H+2cI)u> s.t. dm + m in convex constraint set 
    dm = reshape(convex_sub(model,pm,m,g,Hc),n(1),n(2));
    
    %solve for data augmented ubar, evaluate objective and define gradient
    [objtemp,gtemp,Htemp] = computeGradient(m+dm,model,pm,Q,ssW,Ps,Pr,freq_ind,freq,q,d);
    
    %check for sufficient_decrease before updating m and c
    dec = pm.sigma*(g'*dm(:) + .5*dm(:)'*(Hc.*dm(:)));
    %(objtemp - obj) + abs(dec) %should be < 0
    sd = ( (objtemp - obj) < -abs(dec) );
    if (~sd || dec>0) %dec<=0 should always hold if inner problem solved correctly
        c = min(pm.c2*c,pm.cmax); %increase c (like reducing step size)
        Hc = H + 2*c;
        discarded_oits = discarded_oits + 1;
        disp(['Iteration ' num2str(oit) ' rejected due to insufficient descent, increased c to = ' num2str(c)]);
        fprintf('\n');
    else
        c = max(pm.cmin,c/pm.c1); %decrease c (like increasing step size)
        disp(['decreased c to = ' num2str(c)]);
        obj = objtemp;
        g = gtemp;
        H = Htemp;
        Hc = H + 2*c;
        m_prev = m;
        m = m + dm; %update model
        if (model.redraw == 1)
            ssW = drawWeights(model); %possibly redraw weights
            [obj,g,H] = computeGradient(m,model,pm,Q,ssW,Ps,Pr,freq_ind,freq,q,d);
            Hc = H + 2*c;
        end
        ores = norm(m-m_prev)/(norm(m)+eps); %check how much m is changing
        disp(['ores = ' num2str(ores) ' at outer iteration ' num2str(oit)]);
        fprintf('\n');
    end
    
    %update energy (note that objective changes if ssW changes)
    energy(oit+1) = obj; 
    
    %option to continue iterating (this is mainly for debugging, can remove later)
    if (oit==maxoits+1) %change to (oit==maxoits) to plot intermediate results
        
        %plot objective
        figure(71)
        clf;
        plot(energy); title('objective versus iteration');
        
        %show current velocity estimate
        figure(72)
        clf;
        imagesc(model.xt,model.zt,1./sqrt(m)); 
        title(['current velocity estimate, ' num2str(oit) ' outer iterations']);
        colorbar; caxis([min(model.vmin(:)) max(model.vmax(:))]);
        
        %show update
        figure(62); clf;
        imagesc(model.xt,model.zt,dm);
        colorbar;
        
        %get user input to continue iterating (optional, comment out later)
        mt = input('hit enter to finish or enter number of outer its to continue: ');
        if (~isempty(mt))
            maxoits = maxoits + mt;
        end
    end    
    
    oit = oit + 1;
end
oits = oit-1;
disp([num2str(oits) ' total iterations, ' num2str(oits-discarded_oits) ' accepted']);

end


function W = drawWeights(model)
%redraw weights W for simultaneous sources
W = randn(model.ns,model.nsim);
end


function [obj,g,H] = computeGradient(m,model,pm,Q,ssW,Ps,Pr,freq_ind,freq,q,d)

%notation
N = model.N;
Xint = model.Xint;
Xbnd = model.Xbnd;
mu = pm.mu;
lambda = pm.lambda;

%initialize running sums
obj = 0; %objective
g = zeros(N,1);
H = zeros(N,1);

%loop over frequencies
nf = length(freq);
for v = 1:nf
    w = 2*pi*freq(v);
    %for all sources at once, solve for reduced data-augmented ubar
    src = q(v)*(Ps'*Q*ssW);
    Uv = [sqrt(mu)*lambda*(model.L+model.M(freq(v),m(:)));sqrt(mu)*Pr]\...
        [sqrt(mu)*lambda*src;sqrt(mu)*d{v}*ssW];
    
    %keep running sum of objective obj
    dres = Pr*Uv - d{v}*ssW;
    pres = lambda*(model.L*Uv + model.M(freq(v),m(:))*Uv - src);
    obj = obj + (.5*mu*(dres(:)'*dres(:)) + .5*mu*(pres(:)'*pres(:)))*...
        model.vweights(freq_ind(v));
    
    %keep running sum of gradient g
    Luq = model.L*Uv - src;
    g = g + (sum( real(mu*lambda^2*w*conj(Uv).*(...
        bsxfun(@times,w*Xint,Luq) + ...
        w^3*bsxfun(@times,Xint.*m(:),Uv) - ...
        -1i/2*bsxfun(@times,Xbnd./sqrt(m(:)),Luq) + ...
        w/2*bsxfun(@times,Xbnd.^2,Uv))) ,2))*...
        model.vweights(freq_ind(v));
    
    %keep running sum of diagonal Hessian approximation Hc
    H = H + (sum( mu*lambda^2*w^4*bsxfun(@times,Xint,abs(Uv).^2) + ...
        real(-1i*mu*lambda^2*w/4*bsxfun(@times,Xbnd./(m(:).^(3/2)),Luq).*conj(Uv)) ,2))*...
        model.vweights(freq_ind(v));    
end

end
