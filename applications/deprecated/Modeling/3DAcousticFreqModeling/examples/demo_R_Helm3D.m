%% 3D constant-density acoustic frequency-domain modeling: Basic use of R_Helm3D [LEGACY]
%
% This illustrates the basic use of R_Helm3D. It discretizes the helmholtz 
% operator using the old discretization routine.
% 

%% Setting up the model
% Notice that we need to be careful when choosing PML size and number of points 
% per wavelength. We also need to resize the model (if desired and/or needed) 
% and create the artificial PML layer.
% 
% *Remark:* R_Helm3D requires the velocity to be given in slowness squared.
% 
par.outdir  = '../results/';
par.label   = 'demo_R_Helm3D';
flog = fopen([par.outdir par.label '_log'],'w');
np = parpool_size();

% Define the velocity model
%--------------------------------------
n = [51,51,21];
o = [0 0 0];
d = [50 50 50];
v = ones(n)*(1/(4500^2));

vmin = 1./sqrt(max(v(:)));
vmax = 1./sqrt(min(v(:)));
nlam = 6;
f = 15;

% Computes now the pml width based on max wavelength.
npml = vmax/(f*max(d));
npml = [1 1 1]*ceil(npml);
nt   = n + 2*npml; % Total number of points, with pml included
N    = prod(n);
Nt   = prod(nt);

plog(flog,'* Extending physical domain - creating PML... \n');
Px = opKron(opExtension(n(3),npml(3)),opExtension(n(2),npml(2)),opExtension(n(1),npml(1)));
v = Px*v(:);

plog(flog,'* vmin: ', vmin, ' \n');
plog(flog,'* vmax: ', vmax, ' \n');
plog(flog,'* n:    [ ', n, ' ]  = ', N, ' \n');
plog(flog,'* npml: [ ', npml, ' ], Total: ', Nt, '\n');
plog(flog,'* d:    [ ', d, ' ]\n');
plog(flog,'* f:    ', f, ' Hz \n');
plog(flog,'* nlam: ', nlam, ' \n');

pml.x   = npml(1);
pml.y   = npml(2);
pml.top = npml(3);
pml.bottom = npml(3);

%% Calling R_Helm3D
% 
% Please, check |help R_Helm3D| for more details on the input parameters
% 
plog(flog,'* Discretizing Helmholtz operator with legacy code [FOR COMPARISON PURPOSES ONLY]...  ');
tic; [H,idx]   = R_Helm3D(f,v,ones(length(v),1),o,d,nt,npml,100,2); T = toc;
H = transpose(H);
plog(flog,' done in ', T , ' seconds \n');

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
DScale = 1./sqrt(sum(abs(H).^2,1));
DScale = spdiags(transpose(DScale),0,Nt,Nt);
H = H*DScale; % Normalize matrix

% Create right-hand side
%-------------------------------------
x0 = zeros(Nt,1);
Q  = zeros(Nt,1);
Q(sub2ind(nt, ceil(nt(1)/2) + 1, ...
              ceil(nt(2)/2) + 1, ...
              ceil(nt(3)/2) + 1)) = 1/prod(d);
Q = DScale*Q;


% Distribute Everything
%-------------------------------------
if np > 1
    H  = transpose(H);
    x0 = distributed.zeros(nt);x0 = x0(:);
    Q  = distributed(Q);
    Hr  = distributed(real(H));
    Hi  = distributed(imag(H));
    H   = Hr + 1i*Hi;
    clear Hr Hi;
    spmd,
        dist1 = getCodistributor(x0);
        dist2 = codistributor1d(1,dist1.Partition,[Nt,27]);
        Q     = redistribute(Q,dist1);
        H     = redistribute(H,dist2);
    end
else
    x0 = zeros(Nt,1);
end

% Compute Wavefield
%------------------------------
par.maxit = 2000;
par.tol   = 1e-6;
par.size  = nt;
plog(flog,'\n* Starting Krylov solver\n\n');
if np>1
   tic ; [u,res] = CARPCG(H,idx,Q,x0,par); T = toc;
else
   tic ; [u,res] = CGMN(H,idx,Q,x0,par); T = toc;
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
