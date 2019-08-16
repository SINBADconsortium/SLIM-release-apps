function output = DF(m,Q,input,flag,model,params)
% Frequency domain modeling in the Born approximation. This is the
% Jacobian of F(m,Q,model). Outdated, please use oppDF.m or PDEfunc.m
% instead.
%
% use: 
%   output = DF(m,Q,input,flag,model,{params})
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source and has to be
%                       distributed over the last dimension.
%   input             - flag= 1: vector with gridded slowness perturbation
%                       flag=-1: vectorized data cube of size nrec xnrec x nfreq
%   flag              -  1: forward mode
%                       -1: adjoint mode
%   model.{o,d,n}     - regular grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc
%   model.nb          - number of extra points for absorbing boundary on each side
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array
%
%   params            - optional struct of performance options
%   params.computeLU  - if true, will use LU factors for inverting Helmholtz (default: false)
%   params.nthreads   - if > 0, forces matlab to use multithreading inside spmd blocks (default: 0)
%   params.noext      - if true, does not use the adjoint of opExtension for the adjoint mode, but instead restricts to the computational domain (default: false)
%
%%
% *Example*
% Below example defines a grid [0,1000]^2 with 10 m gridspacing. The
% velocity is 2 km/s. The sources and receivers coincide at x = [0:10:1000]
% and z = 10m. We do the dottest with a randomly generated perturbation.
%
%%
% model.o = [0 0];
% model.d = [10 10];
% model.n = [101 101];        
% model.nb = [10 10];
% model.freq = [5 10 15 20]; 
% model.f0 = 0;
% model.t0 = 0;
% model.zsrc = 10; 
% model.xsrc = 0:10:1000;
% model.zrec = 10; 
% model.xrec = 0:10:1000;
% m = .25*ones(prod(model.n),1);
% Q = speye(length(model.xsrc));
% dm = randn(size(m));
% dD = DF(m,Q,dm,1,model);
% dDt = 0*dD + randn(size(dD));
% dmt = DF(m,Q,dDt,-1,model);
% real(dD'*dDt)
% dm'*dmt
%%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% Updated by: Curt Da Silva, Bas Peters, 2015
% 
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
Nt = prod(nt);
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% Parse performance options
computeLU = false;
nthreads = 0;
extend = true;
if exist('params','var') && ~isempty(params)
    if isfield(params,'computeLU'), computeLU = params.computeLU; end
    if isfield(params,'nthreads'), nthreads = params.nthreads; end
    if isfield(params,'extend'), extend = params.extend; end
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

% distribute frequencies according to standard distribution
freq = model.freq;

% m : s^2/m^2
m = Px*m;

% check source matrix input
if (size(Q,3)==1)&&(isdistributed(Q))
    Q = gather(Q);
end

if nthreads
    set_maxNumCompThreads(nthreads,'automatic');
end

if flag==1
    input = Px*input;
    output = zeros(nsrc*nrec,nfreq);
    % solve Helmholtz for each frequency in parallel
    for k=1:nfreq
        [fm,df] = input2helm_param(m,freq(k)); 
        dfdm = df .* input;
        if size(Q,3)==1, Qk = w(k)*(Ps'*Q); 
        else Qk = w(k)*(Ps'*Q(:,:,k)); 
        end
        % Hk = A * diag( (B*fm).^2 ) + const
        [Hk, A,B] = Helm2D(fm,ot,dt,nt,model.nb);
        
        % Derivative of mapping, fm -> Hk(fm)
        dHk = opMatrix(A)*opDiag((2*B*fm) .* (B*dfdm));            
        
        if computeLU
            [LL UU Pp Qp R]=lu(Hk);
            Uk  = Qp*(UU\(LL\(Pp*(R\(Qk)))));
            dUk  = Qp*(UU\(LL\(Pp*(R\(-dHk*Uk)))));
        else
            % DU[dm] = H[m]\( -dH(m)[dm] * U(m) )   
            Uk = Hk\Qk;                                                 
            dUk = Hk\(-dHk*Uk);
        end                        
        output(:,k) = vec(Pr*dUk);
    end    
else
    output = zeros(prod(model.n),1);
    for k=1:nfreq
        % f : s^2/m^2 -> radians * s/km
        [fm,df] = input2helm_param(m,freq(k)); 
            
        if size(Q,3)==1, Qk = w(k)*(Ps'*Q); 
        else Qk = w(k)*(Ps'*Q(:,:,k)); 
        end
            
        % Hk = A * diag( (B*fm).^2 ) + const
        [Hk, A,B] = Helm2D(fm,ot,dt,nt,model.nb);
            
        % Primary wavefield     
        if computeLU 
            [LL UU Pp Qp R]=lu(Hk);
            Uk       = Qp*(UU\(LL\(Pp*(R\(Qk)))));
            Sk        = -Pr'*reshape(input(:,k),nrec,nsrc);
            V0k       = R'\(Pp'*(LL'\(UU'\(Qp'*Sk))));
        else            
            Uk       = Hk\Qk;      
            Sk        = -Pr'*reshape(input(:,k),nrec,nsrc);             
            V0k       = Hk'\Sk; %Back propagated residual wavefield
        end
        % Gradient of f
        opDf = opDiag(2*(B*fm) .* ( B*df));  
        
        % Adjoint of derivative of mapping, fm -> Hk(fm)
        r = sum( conj(opDf * Uk) .* (A'*V0k) ,2 );
        output = output + real(Pcomp*r);
        end                                         
    end
end
    
  


