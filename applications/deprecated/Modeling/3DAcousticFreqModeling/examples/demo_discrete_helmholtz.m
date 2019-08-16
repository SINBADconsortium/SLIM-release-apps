%% 3D constant-density acoustic frequency-domain modeling: Basic use of discrete_helmholtz
%
% In this example we will present the basic use for the new discretization
% routines, based on the function discrete_helmholtz. This consists of a "black-
% box" function that automatically computes the size of the PML, the smallest
% grid spacing to achieve stability for the chosen frequency, and also deals 
% with interpolation and extrapolation whenever it is needed.
% 
% Naturally, this function internally calls helmholtz_3d.
% 

%% Setting up the model
par.outdir  = '../results/';
par.label   = 'demo_helmholtz_3d';
flog = fopen([par.outdir par.label '_log'],'w');
np = parpool_size();

% Define the velocity model
%--------------------------------------
% The number of points in the model is not important at all!
% It will be resized to the smallest size satisfying the stability
% condition inside discrete_helmholtz.
% What truly matters is the physical size (that is, nv*dv)
% 
model.nv   = [51,51,21];
model.ov   = [0 0 0];
model.dv   = [50 50 50];
model.v    = ones(model.nv)*4500;
model.unit = 'm/s';


%% Calling discrete_helmholtz
% This is all needed: the model structure and the target frequency.
% An optional structure |opts| can be passed, allowing the user to impose
% some conditions (as for instance, imposing a certain number of points in the 
% PML, or use a different number of points per wavelength) but we don't use it
% in this example. Please, check the documentation by running 
% |help discrete_helmholtz| inside Matlab to a full description on this 
% structure.
% The last parameter, "flog" is just a file log. It is as well optional. You 
% can either ignore it or pass 0 - no messages will be printed in this case.
helm = discrete_helmholtz(model,15,[],flog);

%% Computing wavefield
% We already have the discrete Helmholtz operator; now we just need to define
% the right hand side and then call a solver.
% 
% Since the operator is stored in band-storage format, we will use CARPCG to
% solve the wave-equation. If you wish to use another solver instead, you will
% need to convert the operator from band-storage format to Matlab's sparse 
% matrix format. This is achieved with
% 
% |A = H2sparse(helm.coef,helm.idx);|
% 
% and then using the desired solver normally.

% Normalize rows of the matrix - This is needed only because we will use CARPCG
% to solve the linear system; Otherwise this can/should be skipped
Nt = prod(helm.nt);
DScale = 1./sqrt(sum(abs(helm.coef).^2,1));
DScale = spdiags(transpose(DScale),0,Nt,Nt);
helm.coef = helm.coef*DScale; % This is actually a multiplication from the LEFT!

% Create right-hand side (point source)
%---------------------------------------
x0 = zeros(Nt,1);
Q  = zeros(Nt,1);
Q(sub2ind(helm.nt, ceil(helm.nt(1)/2) + 1,...
                   ceil(helm.nt(2)/2) + 1,...
                   ceil(helm.nt(3)/2) + 1)) = 1/prod(helm.d);
Q = DScale*Q;

% Distribute Everything
%-------------------------------------
if np > 1
   plog(flog,'\n* Distributing operator over ', np, ' workers...');
   H  = transpose(helm.coef);
   x0 = distributed.zeros(helm.nt);x0 = x0(:);
   Q  = distributed(Q);
   Hi  = distributed(imag(H));
   H   = distributed(real(H));
   H   = H + 1i*Hi;
   spmd,
      dist1 = getCodistributor(x0);
      dist2 = codistributor1d(1,dist1.Partition,[Nt,length(helm.idx)]);
      Q     = redistribute(Q,dist1);
      H     = redistribute(H,dist2);
   end
   helm.coef = transpose(H);
   clear H Hi;
else
    x0 = zeros(Nt,1);
end

% Compute Wavefield
%------------------------------
par.maxit = 2000;
par.tol   = 1e-5;
par.size  = helm.nt;
plog(flog,'\n* Starting Krylov solver\n\n');
if np>1
   tic ; [u,res] = CARPCG(helm.coef,helm.idx,Q,x0,par); T = toc;
else
   tic ; [u,res] = CGMN(helm.coef,helm.idx,Q,x0,par); T = toc;
end
plog(flog,'* Krylov stopped after ', T , ' seconds \n');

% Write Results
%------------------------------
outfile1 = [par.outdir par.label '_out.mat'];

if np>1
    u = gather(u);
    Q = gather(Q);
end

plog(flog,'* All done... saving results... \n');
plog(flog,'* Writing ', outfile1, '... \n');
save(outfile1,'res','Q','u');

fclose(flog);

