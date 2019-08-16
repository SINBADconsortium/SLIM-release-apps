addpath(genpath('../'));
model.o = [0 0];
model.d = [10 10];
model.n = [101 101];
model.nb = [20 20];
model.freq = [10 15];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 15;
model.xsrc = 0:10:1000;
model.zrec = 10;
model.xrec = 0:5:1000;
v0 = 2000;
model.m = 1e6/v0.^2*ones(prod(model.n),1);

Q = speye(length(model.xsrc));
m = model.m;
params = struct; 
params.computeLU = true;
params.extend    = 1;
[D0, J0] = F(m,Q,model,params);
D0 = gather(D0);

%Randomize over perturbations
nTrials = 1;

%Error in order of convergence tolerance, since the order of convergence is asymptotic + limited precision from having to solve linear systems
tol = 0.15;

% if h << 1e-6, you run in to floating point issues so the test might fail then
h = 10.^(-6:0);

for j=1:nTrials

    dm = randn(model.n); 
%     dm(1,:) = 0; dm(end,:) = 0;
%     dm(:,1) = 0; dm(:,end) = 0; 
    dm = dm(:);
    
    dD = gather( J0*dm );       
    
    % Gradient test    
    e0 = zeros(length(h),1);
    e1 = zeros(length(h),1);
    for i=1:length(h)
        D1 = gather(F(m+h(i)*dm,Q,model,params));
        e0(i) = norm(D1-D0,'fro');
        e1(i) = norm(D1-D0-h(i)*dD,'fro');
    end
    
    % Order of convergence
    h0 = log10(mean(e0(2:end)./e0(1:end-1)));
    h1 = log10(mean(e1(2:end)./e1(1:end-1)));
    
    % First order convergence without gradient
%    assert( abs( h0 - 1 ) < tol );
    
    % Second order convergence with gradient
%    assert( abs( h1 - 2 ) < tol );
end

figure;loglog(h,e0);hold on;loglog(h,h);loglog(h,h.^2);loglog(h,e1);legend('e0','h','h2','e1')