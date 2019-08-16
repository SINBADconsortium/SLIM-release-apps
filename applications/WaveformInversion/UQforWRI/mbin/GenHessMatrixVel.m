function [V G AA EIG Px] = GenHessMatrixVel(m,Q,D,model,params)
% GenHessMatrix - A function that compute the approximate Hessian for WRI around 
% the optimal solution
% 
% 
% Usage:
%   [V G EIG] = GenHessMatrixVel(m,Q,D,model,params)
%
% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Input:
%   m                 - vector with gridded velocity in [m/s]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source
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
%   params.extend     - if true, will use true adjoint of extension operator to %                       return to the computational grid (default: false)
%
% Output:
%   [V G EIG]
%   

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
% z = zeros(size(Pr,1),1);
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

%if ~isempty(dm), dm = Px*dm; end

if nthreads
    set_maxNumCompThreads(nthreads,'automatic');
end

if isfield(params,'WS')
    nsrc = size(params.WS,2);
else
    nsrc = size(Q,2);
end

G    = zeros(length(m), nsrc, nfreq);
V    = zeros(size(Prsp,1), length(m), nfreq);
EIG  = zeros(size(Prsp,1),nfreq);
for k=1:length(freq)
    if length(lambdas)>1
        lambda = lambdas(k);
    else
        lambda = lambdas;
    end
    [fm,dfm,ddfm] = input2helm_param_vel(m,freq(k));
    if size(Q,3)==1, Qk = w(k)*(Ps'*Q);
    else Qk = w(k)*(Ps'*Q(:,:,k));
    end
    
    % Hk = A * diag( (B*fm).^2 ) + const
    [Hk, A1,B1] = Helm2D(fm,ot,dt,nt,model.nb);
%    M = (B1.*B1) * (2*pi*freq(k)*1e-3)^2;
    M = opDiag_swp(2*(B1*fm) .* ( B1*dfm ));
    
    % set up and solve the data-augmented wave equation
    A       = [lambda*Hk;Prsp];
    if isfield(params,'WR')
        S       = [lambda*Qk;params.WR*D(:,:,k)];
    else
        S       = [lambda*Qk;D(:,:,k)];
    end
    
    if isfield(params,'WS')
        S = S * params.WS;
    end
    
    GG      = A'*A;
   
    U        = GG\(A'*S);
    G(:,:,k) = M * U;
    AA{k}    = A1;
    
    U2                    = full((Hk' \ Prsp')');
    [UU S Vloc]           = svd(U2,'econ');
    V(:,:,k)              = Vloc';
    EIG(:,k)              = diag(S);
    
   
end


end
