function [ varargout ] = PDEfunc( func, m,Q,input,Dobs,model,params )
%PDEFUNC An abstract function that computes some output quantity depending
% on a baseline model, a source matrix, and possibly some input
% perturbation + data. Used by misfit_setup, F.m, oppDF, oppHGN, oppH. 
%
% Usage:
%   varargout = PDEfunc( func, m, Q, {input}, {Dobs}, model, params);
%
% Author: Curt Da Silva, Bas Peters, Tristan van Leeuwen 
%
% Input:
%   func               - string, one of
%                         - PDEopts.FORW_MODEL    - forward modelling operator
%                         - PDEopts.JACOB_FORW    - Jacobian of PDEopts.FORW_MODEL
%                         - PDEopts.JACOB_ADJ     - Jacobian adjoint of PDEopts.FORW_MODEL
%                         - PDEopts.HESS_GN       - Gauss-Newton hessian of PDEopts.FORW_MODEL
%                         - PDEopts.HESS          - Hessian of 'forw_model'
%                         - PDEopts.OBJ           - least-squares objective + gradient
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source
%   input             - input model perturbation, used by 'jacob_forw','hess_gn','hess'
%   Dobs              - observed data, used by 'hess','obj'
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%
%   params            - optional struct of performance options
%   params.pde_solve  - string for which method to solve the Helmholtz equation, one of 
%                         - 'iter' - (default) uses iterative methods 
%                                  .ls_tol    - relative residual tolerance (default: 1e-6)
%                                  .ls_solver - one of LinSolveOpts.SOLVE_CGMN, LinSolveOpts.SOLVE_CRMN (default), LinSolveOpts.SOLVE_FGMRES
%                                  .ls_niter  - max number of iterations (default: 2000)
%                        
%                         - 'lu'         - uses Matlab's built in sparse LU decomposition
%                         - 'backslash'  - uses Matlab's backslash, not recommended
%   params.nthreads   - if > 0, forces matlab to use multithreading inside spmd blocks (default: 0)
%   params.extend     - if true, will use true adjoint of extension operator to return to the computation grid (default: false) 
%
%
%
% Output:
%   if func=='obj', 
%      {obj, gradient}
%   
%   if func=='forw_model','jacob_forw'
%      forward, born-scattering wavefield, resp., of size [nrec x nsrc x nfreq]
%
%   if func=='jacob_adj','hess_gn','hess'
%      migrated image, gauss-newton hessian, hessian image, resp., of size prod(model.n) x 1
%
%

% comp. grid
b =1;
FORW_MODEL = PDEopts.FORW_MODEL; JACOB_FORW = PDEopts.JACOB_FORW; JACOB_ADJ = PDEopts.JACOB_ADJ;
HESS_GN = PDEopts.HESS_GN; HESS = PDEopts.HESS; OBJ = PDEopts.OBJ; OP1 = PDEopts.OP1;
modes = {OBJ,FORW_MODEL,JACOB_FORW,JACOB_ADJ,HESS_GN,HESS,OP1};

if isempty(find(ismember(func,modes))), error('Must choose a preset function.'); end

ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% Parse performance options
SOLVE_ITER = 'iter'; SOLVE_LU = 'lu'; SOLVE_BACKSLASH = 'backslash';
pde_solve = SOLVE_LU;
nthreads = 0;
extend = false;

if exist('params','var') && ~isempty(params)    
    if isfield(params,'pde_solve'), pde_solve = params.pde_solve; end
    if isfield(params,'nthreads'), nthreads = params.nthreads; end
    if isfield(params,'extend'), extend = params.extend; end    
    % Default iterative solver options
    if strcmp(pde_solve,SOLVE_ITER) 
        lsopts = LinSolveOpts();
        if isfield(params,'ls_tol'),lsopts.tol = params.ls_tol; else lsopts.tol = 1e-6; end
        if isfield(params,'ls_solver'),lsopts.solver = params.ls_solver; else 
            lsopts.solver = LinSolveOpts.SOLVE_CRMN; end
        if isfield(params,'ls_iter'), lsopts.maxit = params.ls_iter; else lsopts.maxit = 2000; end
        if isfield(params,'precond'), lsopts.precond = params.precond; end   
    end    
    if isfield(params,'helm_scheme'), helm_scheme = params.helm_scheme; else helm_scheme = PDEopts.HELM2D_JO9P; end
else
    params = [];
end

% define wavelet
w = fwi_wavelet(model.freq,model.t0,model.f0);

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
Pe = opKron(opExtension(model.n(2),model.nb(2),0),opExtension(model.n(1),model.nb(1),0));

if extend, Pcomp = Px'; else Pcomp = Pe'; end

freq = model.freq;

% m : s^2/km^2
mx = Px*m;
if ~isempty(input),
    if strcmp(func,JACOB_ADJ)
        input = reshape(input,nsrc*nrec,nfreq);
    else
        input = Px*input; 
    end
end
if ~isempty(Dobs)
    Dobs = reshape(Dobs,nrec,nsrc,nfreq); 
end
% check source matrix input
if (isdistributed(Q))
    error('PDEfunc can only handle non-distributed inputs');
end

if nthreads
    set_maxNumCompThreads(nthreads,'automatic');
end

% Allocate output
switch func 
    case OBJ
        f = 0; g = zeros(prod(model.n),1);
    case {FORW_MODEL,JACOB_FORW}
        output = zeros(nsrc*nrec,nfreq);
    otherwise
        output = zeros(prod(model.n),1);
end

for k=1:nfreq
     % f : s^2/(km)^2 -> radians * s/m
     
     [fm,df,ddf] = input2helm_param(mx,freq(k));
      
     % Hk = A*diag((B*fm).^2) + const
     switch helm_scheme
       case PDEopts.HELM2D_JO9P        
         [Hk,A,B] = Helm2D(fm,ot,dt,nt,model.nb);
       case PDEopts.HELM2D_CHEN9P
         [Hk,A,B] = Helm2D_opt(fm,dt,nt,model.nb,freq(k),model.f0);
       otherwise
         error('Unrecognized helmholtz scheme');
     end
     
     % Per-frequency source
     if size(Q,3)==1, Qk = w(k)*(Ps'*Q);
     else Qk = w(k)*(Ps'*Q(:,:,k)); end
     
     % dM = d/dm diag((B*fm).^2);
     dM = opDiag(2*(B*fm) .* ( B*df ));
     
     switch func
         case {JACOB_FORW,HESS_GN,HESS}
            dHk = opMatrix(A) * dM * opDiag(input);                    
     end     
     
     switch pde_solve
         case SOLVE_ITER
             [R,idx] = mat2R(Hk); R = R.'; idx = idx.';                           
             helm_opts = struct; helm_opts.mode = opBandStorage.MULT_DIV; 
             helm_opts.solve_opts = lsopts;
             Hk = opBandStorage(R,idx,helm_opts);
             Hinv = @(x) Hk\x; Hinvadj = @(x) Hk'\x;
         case SOLVE_LU             
             [LL,UU,Pp,Qp,R] = lu(Hk);
             Hinv = @(x) Qp*(UU\(LL\(Pp*(R\(x))))); 
             Hinvadj = @(x) R'\(Pp'*(LL'\(UU'\(Qp'*x))));
         case SOLVE_BACKSLASH
             Hinv = @(x) Hk\x;
             Hinvadj = @(x) Hk'\x;     
     end     
     
     % Every method uses the background wavefield
     Uk  = Hinv(Qk);
     
     switch func,  
         % Least squares objective + gradient
         case OBJ
             SrcEst = 0;
             if isfield(params, 'SrcEst')
                SrcEst = params.SrcEst;
             end
             
             
             
             if SrcEst > 0
                dtmp1         = vec(Pr*Uk);
                dtmp2         = vec(Dobs(:,:,k));
                src_alpha  = dtmp1'*dtmp2 / norm(dtmp1)^2;
                Uk                = src_alpha * Uk; 
             end
             
             r     = Pr*Uk - squeeze(Dobs(:,:,k));
             f     = f + 0.5*norm(r,'fro')^2;
             if nargout >= 2
                Vk    = Hinvadj( -Pr'*r );
                g     = g + real(Pcomp*sum( conj(dM * Uk) .* (A'*Vk) ,2 ) );
             end
             
         % Output is just the forward wavefield at the receivers
         case FORW_MODEL
             output(:,k) = vec(Pr*Uk);    
             
         % Output is just the derivative wavefield at the receivers
         case JACOB_FORW
             dUk = Hinv(-dHk*Uk);
             output(:,k) = vec(Pr*dUk);
             
         % Output is just the migrated image
         case JACOB_ADJ        
             Vk = Hinvadj(-Pr'*reshape(Dobs(:,:,k),nrec,nsrc));
             r = sum( conj(dM * Uk) .* (A'*Vk) ,2 );
             output = output + real(Pcomp*r);
       case OP1             
             Vk = Hinvadj(-Pr'*squeeze(Dobs(:,:,k)));
             U1 = sum(Uk,2)*(2*pi*freq(k))^(-1); V1 = sum(Vk,2)*(2*pi*freq(k))^(-1);
             dz = opFunction(nt(1),nt(1),@(x,mode) stencil_mvp1d(x,{@(x) -1/(2*dt(1))*x,@(x) 0*x,@(x)1/(2*dt(1))*x},mode,nt(1)));
             dx = opFunction(nt(2),nt(2),@(x,mode) stencil_mvp1d(x,{@(x) -1/(2*dt(2))*x,@(x) 0*x,@(x)1/(2*dt(2))*x},mode,nt(2)));
             Dz = opKron(opDirac(nt(2)),dz); Dx = opKron(dx,opDirac(nt(1)));             
             out = (Dz*U1).*conj(Dz*V1) + (Dx*U1).*conj(Dx*V1);
             output = output + real(Pcomp*out);  
             
         % The GN Hessian is the composition of 'jacob_forw' and 'jacob_adj'    
         case HESS_GN
             dUk = Hinv(-dHk*Uk);
             dUk = Hinvadj(-Pr'*Pr*dUk);
             r = sum( conj(dM * Uk) .* (A'*dUk) ,2 );
             output = output + real(Pcomp*r);
             
         % Original gradient is conj(V1*V2).*(A'*V3), so the Hessian is the
         % derivative of this quantity
         case HESS
             dUk = Hinv(-dHk*Uk);
             V1 = dM;     dV1 = opDiag_swp((2*(B*df).^2 .* input + (2*B*fm) .* (B * (ddf.*input))));             
             V2 = Uk;     dV2 = dUk;
             V3  = -Hinvadj( Pr'*Pr * Uk - Pr'*Dobs(:,:,k) );
             dV3 = -Hinvadj( Pr'*Pr * dUk + dHk' * V3  );       
             r = real(sum(conj(dV1*V2 + V1*dV2) .* (A'*V3) + conj(V1*V2) .*(A'*dV3) ,2));            
             output = output + real(Pcomp*r);
     end    
end
if strcmp(func,OBJ)
    varargout = {f, g};
else
    varargout{1} = output;
end

end
