%% 3D constant-density acoustic frequency-domain modeling: Testing Helmholtz derivative
%
% In this test we certify that the function helmholtz_3d_derivative in fact 
% computes the numerical derivative of the function helmholtz_3d. Notice that 
% this is crucial for computing the action of the Jacobian (without explicitly 
% forming the Jacobian).
%
% *Author:* <mailto:rlago@eos.ubc.ca Rafael Lago>

%% Theory
% From Taylor's theorem we know that
% 
% $H(m) = H(m+h\delta m) + h\partial H(m+h\delta m)\delta m + \mathcal O(h^2)$
% 
% That means that the second order perturbation
% 
% $H(m) - H(m+h\delta m) - h\partial H(m+h\delta m)\delta m$ (1)
% 
% should be of the same order of magnitude as $\mathcal O(h^2)$ for any model 
% $m$ and a small enough perturbation $\delta m$. If the model is given in 
% slowness squared rather than meters per second, the second derivative should
% be zero and thus the first order approximation should be exact up to machine
% precision. In this test we concern the meters per second unit only, but with
% little modification, this script also shows that helmholtz_3d_derivative
% also represents the derivative of helmholtz_3d when using slowness squared.
% 
% In this test, we will consider the simple 3D Camembert velocity model as $m$
% and a random perturbation $\delta m$. Instead of measuring $H(m)-H(m+ \delta m)$ 
% however, we will sample it with a vector $p$, as well chosen at random, and 
% the test we effectively perform is
% 
% $H(m)p - H(m+h\delta m)p - h(\partial H(m+h\delta m)\delta m)p$ (2)
% 
% for varying values of $h$, from $10^{-1}$ to $10^{-10}$. We also impose that 
% $||\delta m||=||m||$ where $||.||$ denotes the Euclidean norm of the 
% vectorized model $m$.
% 
%

%% Computational details
% 
% The actual derivative of helmholtz_3d should be a tensor of dimension 
% $(n_x\times n_y\times n_z)^3$ and therefore the storage of such a
% structure is considerably expensive. However, it is possible to compute and
% store the tensor dH/dm using band-storage format. The total storage cost for
% the tensor $\partial H/\partial m$ computed by helmholtz_3d_derivative is 
% exactly the same as that of the matrix H, computed using helmholtz_3d.
% 
% We refer to the technical report [1] for more details on how to implement the
% tensor and the tensor-vector product using the band-storage format.
% 
% *Remark:* This is specifically the tensor for the 27 points stencil. The 
% implementation of the tensor for the simpler 7 points stencil will be provided
% in a future release.
%


%% Test case
% For this test we use the simple 3D camembert model. It consists of
% a homogeneous model with a high velocity ellipsoid in the center.
% The following figure shows a slice of the model.

clear
flog = fopen(['../results/test_helmholtz_derivative_log'],'w');


% Read velocity model (and plot a slice)
%--------------------------------------------------------------------
vfile      = [ '../data/m_true.rsf' ];
[model.v model.nv model.dv model.ov]  = rsf_read_all(vfile);
model.unit = 'm/s';
plot_slice(model.v,11,'z');


%% Script
% Follows the test script. It prints some details of the velocity model being
% used and the progress of the numerical test for each value of h.

% Obtains a stable computational grid 
% (does not discretize yet)
%--------------------------------------------------------------------
opts.discretize=false;
cg = discrete_helmholtz(model,6,opts,flog);

N = prod(cg.nt);

% Choose the random perturbation dm and a random vector p
%---------------------------------------------------------
h  = 10.^([-1:-1:-10]);
rng(123456);
dm = randn(model.nv);
dm =(dm/norm(dm(:))) *norm(model.v(:));
m  = cg.pg2cg(model.v);
dm = cg.pg2cg(dm);
p  = randn(N,1);

% Discretize for H(m)
%---------------------
[H0,idx] = helmholtz_3d(m,cg.d,cg.pml,cg.f,model.unit);
H0p = Hmvp(H0,idx,p);

plog(flog,'Testing the partial derivative of the Helmholtz operator w.r.t. m\n');

% Loop over h
%-------------
for i=1:length(h)
   [H, idx] = helmholtz_3d(m + h(i)*dm,cg.d,cg.pml,cg.f,model.unit);
   Hp   = Hmvp(H,idx,p);
   
   [dH, ~ ] = helmholtz_3d_derivative(m + h(i)*dm,cg.d,cg.f,model.unit);
   dH       = H2sparse(dH,idx);
   dmd      = spdiags(dm(:),0,N,N);

   diff1(i) = norm(Hp - H0p);
   diff2(i) = norm(Hp - H0p - h(i)*(dH*dmd)*p);
   plog(flog,'Order of perturbation   - ||h*dm||:                               ',h(i),'\n');
   plog(flog,'First order error       - ||H(m) - H(m+h*dm)||:                   ',diff1(i),'\n');
   plog(flog,'Second order correction - ||h*dH(m+h*dm)*dm||:                    ',norm(h(i)*(dH*dmd)*p),'\n');
   plog(flog,'Second order error      - ||H(m) - H(m+h*dm) - h*dH(m+h*dm)*dm||: ',diff2(i),'\n\n');
end


%% Results
% As previously mentioned, we expect the second order error (1) to 
% behave like $h^2$. The following plot clearly shows that the as h decreases,
% the second order approximation follows the same behaviour (i.e. the lines are
% parallel) until reaching machine precision, which is $2.2\times 10^{-16}$.
% Notice that any further computation is not actually meaningful, and we observe
% that the second order error turns into a flat line after that point. This 
% behaviour is expected.
% 
% Likewise, we show the first order error as well. As expected, it behaves like 
% $h$ (i.e. the lines are parallel).
% 

figure;
loglog(h,diff1,h,h,'k--');
xlabel('h');ylabel('error');
hold on;
loglog(h,h.^2,'g--');
hold on;
loglog(h,diff2,'r');
xlim([h(end) h(1)]);
legend('First order error','h','h^2','Second order error','Location','NorthWest');

fclose(flog);

%% References
% <2014-Lago-Efficient_computation.pdf [1]> Lago, R. 2014 - Efficient 
% computation of dH/dm tensor for a 27-points stencil with mass lumping
%
