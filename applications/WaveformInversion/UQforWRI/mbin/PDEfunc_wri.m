function varargout = PDEfunc_wri(func,m,Q,dm,D,model,params)
% HGN_PEN - An abstract function that computes some output quantity
% depending on a baseline model, source matrix, and possibly some input
% perturbation + data. Used by misfit_setup, misfit_pen.m, oppHGNpen, oppHpen
% 
% Usage:
%   varargout = PDEfunc( func, m, Q, {input}, D, model, params );
%
% Author: Curt Da Silva, Bas Peters, Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Input:
%   func              - string, one of
%                        - 'obj'          - WRI least squares objective
%                        - 'hess_gn'      - Gauss-Newton hessian of 'obj'
%                        - 'hess'         - Hessian of 'obj'
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source
%   input             - input perturbation, used by 'hess_gn','hess'
%   D                 - observed data
%
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%  
%   params            - struct of parameters
%   params.lambda     - penalty parameter for WRI
%   params.extend     - if true, will use true adjoint of extension operator to %                       return to the computational grid (default: true)
%
% Output:
%   [f,g,h,w,f_aux]   - if func == 'obj'
%   
%   otherwise, Hessian, Gauss-Newton Hessian applied to the input perturbation
%
warning off;
OBJ = 'obj'; HESS_GN = 'hess_gn'; HESS = 'hess';
modes = {OBJ,HESS_GN,HESS};

if isempty(find(ismember(func,modes))), error('Must choose a preset function.'); end

ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
Nt = prod(nt);
[zt,xt] = odn2grid(ot,dt,nt);

lambdas = params.lambda;
assert(length(lambdas)==1 || length(lambdas)==length(model.freq),'Need a single lambda or one per frequency');


% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% Parse performance options
nthreads = 0;
extend = true;
C = ones(nrec,nsrc); 
if exist('params','var') && ~isempty(params)
    if isfield(params,'nthreads'), nthreads = params.nthreads; end
    if isfield(params,'extend'), extend = params.extend; end
    if isfield(params,'C'), C = params.C; end
end

% define wavelet
w = fwi_wavelet(model.freq,model.t0,model.f0);

freq = model.freq;

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
Pe = opKron(opExtension(model.n(2),model.nb(2),0),opExtension(model.n(1),model.nb(1),0));

if extend, Pcomp = Px'; else Pcomp = Pe'; end
idx_m = logical( Pe*(ones(size(Pe,2),1)) );

% check source matrix input
if (size(Q,3)==1)&&(isdistributed(Q))
    Q = gather(Q);
end

% Get the restriction operator in sparse matrix format
e = zeros(size(Pr,1),1); e(1) = 1;
I = []; J = []; V = [];
z = zeros(size(Pr,1),1);
% for i=1:size(Pr,1)
%     if i > 1, e(i-1) = 0; end, e(i) = 1;
%     z = Pr'*e;
%     J1 = find(z);
%     I = [I;i*ones(length(J1),1)]; J = [J;vec(J1)]; V = [V;vec(z(J1))];
% end
% Prsp = sparse(I,J,V,size(Pr,1),size(Pr,2));
Prsp   = sparsedouble(Pr);
if isfield(params,'WR')
    Prsp = params.WR * Prsp;
end

D = reshape(D,[nrec,nsrc,nfreq]);
m = Px*m;

if ~isempty(dm), dm = Px*dm; end

if nthreads
    set_maxNumCompThreads(nthreads,'automatic');
end

% Allocate output
if strcmp(func,OBJ)
    f = 0; g = 0; h = [];
    f_aux = struct; f_aux.pde = 0; f_aux.dat = 0;
else
    output = zeros(prod(model.n),1);
end
A1_2 = [];
for k=1:length(freq)
%labindex
    if length(lambdas)>1
        lambda = lambdas(k);
    else
        lambda = lambdas;
    end
    
    
    if isfield(model,'sigma');
        sigma_k = model.sigma(k);
    end
    [fm,dfm,ddfm] = input2helm_param(m,freq(k));
    if size(Q,3)==1, Qk = w(k)*(Ps'*Q);
    else Qk = w(k)*(Ps'*Q(:,:,k));
    end
    
    % Hk = A * diag( (B*fm).^2 ) + const
    [Hk, A1,B1] = Helm2D(fm,ot,dt,nt,model.nb);
    B2 = (B1.*B1) * (2*pi*freq(k)*1e-3)^2;
    
    % set up and solve the data-augmented wave equation
    if isfield(params, 'SIGMA_h')
        A       = [lambda*Hk; params.SIGMA_h*Prsp];
    else
        A       = [lambda*Hk;Prsp];
    end
    % S       = [lambda*Qk;D(:,:,k)];
    if isfield(params,'WR')
        if isfield(params, 'SIGMA_h')
            S       = [lambda*Qk;params.SIGMA_h * params.WR*D(:,:,k)];
        else
            S       = [lambda*Qk;params.WR*D(:,:,k)];
        end
    else
        if isfield(params, 'SIGMA_h')
            S       = [lambda*Qk;params.SIGMA_h * D(:,:,k)];
        else
            S       = [lambda*Qk;D(:,:,k)];
        end
            
    end
    
    if isfield(params,'WS')
        S = S * params.WS;
    end
    
    %% For Gibbs Sampling
    if isfield(params,'MCMC')
        if params.MCMC > 0
            S = S + randn(size(S));
        end
    end
    
    
    G       = A'*A;
    
    U       = G\(A'*S);
    Dp      = Prsp*U;
    if isfield(params,'WS')
        pde_res = Hk*U-Qk*params.WS;
    else
        pde_res = Hk*U-Qk;
    end
    
    % if isfield(params, 'MCMC')
    %     if params.MCMC > 0
    %         pde_res = pde_res - randn(size(pde_res)) / lambda;
    %     end
    % end
        
    
    M = opDiag_swp(2*(B1*fm) .* ( B1*dfm ));
    
    switch func
        case OBJ
            
            Dobsk = D(:,:,k);
            if isfield(params,'WR')
                Dobsk = params.WR * Dobsk;
            end
            if isfield(params,'WS')
                Dobsk = Dobsk * params.WS;
            end
        
            fdata = .5*norm(C.*(Dp - Dobsk),'fro')^2;
            
            
            fpde  = .5*lambda^2*norm(pde_res,'fro')^2;
            if isfield(model,'sigma')
                fdata = fdata / sigma_k^2;
                fpde  = fpde / sigma_k^2;
            end
            
            f     = f + fdata + fpde; f_aux.pde = f_aux.pde + fpde; f_aux.dat = f_aux.dat + fdata;
            V1    = M * conj(U);      V2    = A1'*pde_res;      
            if isfield(model,'sigma')
                g     = g + lambda^2 * Pcomp * (sum(real(B1'*(V1.*V2)),2)) / sigma_k^2;
            else
                g     = g + lambda^2 * Pcomp * (sum(real(B1'*(V1.*V2)),2));
            end
            
            
            if isempty(A1_2), A1_2 = A1'*A1; %A1_2 = speye(size(A1_2,2));
                h = sparse([],[],[],prod(model.n),prod(model.n),nnz(A1_2)); 
            end
            
            if isfield(params,'WS')
                nsrc = size(params.WS,2);
            else
                nsrc = size(Q,2);
            end
            
            for i=1:nsrc
               P = M*U(:,i); P = spdiags(P(idx_m),0,prod(model.n),prod(model.n));
               if isfield(model,'sigma')
                   h = h + lambda^2 * real( P'*A1_2(idx_m,idx_m) * P ) / sigma_k^2;
               else
                   h = h + lambda^2 * real( P'*A1_2(idx_m,idx_m) * P );
               end
            end

        case HESS_GN
            dHdm = opMatrix(A1)*M*opDiag_swp(dm);
            dU = lambda^2 * ( G \ (-dHdm'*(pde_res) - Hk' * (dHdm*U)));
            Z1 = Prsp*dU; Z2 = lambda * (dHdm * U + Hk * dU);
            Z1 = Prsp' * Z1 + lambda * Hk' * Z2;
            Z1 = lambda^2 * ( G \ Z1 );
            R = conj( -M'*A1' * (pde_res)) .* Z1;
            R = sum(R - conj(M*U)  .* (A1' * Hk * Z1),2);
            Z2 = sum( conj( M * U ) .* (A1' * (lambda*Z2)) , 2);
            if isfield(model,'sigma')
                output = output + Pcomp*(real(R + Z2)) / sigma_k^2;
            else
                output = output + Pcomp*(real(R + Z2));
            end
            
        case HESS
            dHdm  = opMatrix(A1)*M*opDiag_swp(dm);
            dM    = opDiag_swp((2*(B1*dfm).^2 .*dm + (2*B1*fm) .* (B1 * (ddfm).*dm)));
            dU    = lambda^2 * (G \ (-dHdm'*(pde_res) - Hk' * (dHdm*U)));
            V1    = M*conj(U);         dV1 = dM*conj(U) + M*conj(dU);
            V2    = A1'*pde_res;       dV2 = A1'*( dHdm * U + Hk*dU);
            if isfield(model,'sigma')
                output = output + lambda^2*Pcomp*(sum(real(B1'*( dV1.*V2 + V1.*dV2)),2)) / sigma_k^2;
            else
                output = output + lambda^2*Pcomp*(sum(real(B1'*( dV1.*V2 + V1.*dV2)),2));
            end
    end
end

if strcmp(func,OBJ)
    varargout = {f,g,h,f_aux};
else
    varargout{1} = output;
end

end
