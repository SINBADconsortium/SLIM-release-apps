function [f,g,w] = JI_src(m,Q,D,I,model,params)
% misfit for waveform inversion with source estimation
%
% use:
%  [f,g,w] = JI_src(m,Q,D,I,model,params)
%
% input:
%   m - vector with gridded squared-slowness [km^2/s^2] (see also F)
%   Q - matrix containing gridded source functions (see also F)
%   D - observed data as distributed vector (consistent with output of F(m,Q,model))
%   I - index set indicating which sources to evaulate the misfit on
%   model - struct with model parameters (see also F)
%   params.misfit  - function handle to penalty function of the form 
%                    [f,g,h] = misfit(x) where f is a scalar, g is the
%                    gradient of the misfit w.r.t. x and h is (an
%                    approximation) of the hessian w.r.t. x (as vector).
%                    default: @(x)twonorms(x).
%   params.src_est - function handle to misfit used for source estimation
%                    default: params.misfit. Set params.src_est = [] for no
%                    source estimation.
%   params.C       - source-receiver weighting matrix of size nrec x nsrc.
%                    default: no weighting.
%
% output:
%   f - misfit value
%   g - gradient of misfit w.r.t. m
%   w - estimated source wavelet.

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if nargin < 6
    params = [];
end

% dimensions of data
nfreq = length(model.freq);
nrec  = length(model.xrec)*length(model.zrec);
nsrc  = size(Q,2);

% default values
misfit  = @(x)twonorms(x);
src_est = @(x)twonorms(x);
C       = ones(nrec,nsrc);

% check options
if isfield(params,'misfit')
    misfit = params.misfit;
end
if isfield(params,'src_est')
    src_est = params.src_est;
end
if isfield(params,'C')
    C = params.C;
end

% select sources, redefine Q,C,nsrc accordingly
%P  = oppBlockDiag(nfreq,opKron(opRestriction(nsrc,I),opDirac(nrec)));
P = oppKron2Lo(opKron(opDirac(nfreq),opRestriction(nsrc,I)),opDirac(nrec));
D  = P*D;
Q  = Q(:,I);
C  = C(:,I);
nsrc = length(I);

% model data
[Dt,Jt] = F(m,Q,model);

% define mask and apply
CC = oppBlockDiag(nfreq,opDiag(C(:)));
D = pSPOT.utils.distVectorize(data_redistribute(D,nrec,nsrc,nfreq,'freq'));
Dt = pSPOT.utils.distVectorize(data_redistribute(Dt,nrec,nsrc,nfreq,'freq'));
D  = CC*D;
Dt = CC*Dt;
D = pSPOT.utils.distVectorize(data_redistribute(D,nrec,nsrc,nfreq,'srcfreq'));
Dt = pSPOT.utils.distVectorize(data_redistribute(Dt,nrec,nsrc,nfreq,'srcfreq'));


% estimate source wavelet and define weight for data
w = ones(nfreq,1);
if ~isempty(src_est)
    fw = @(w)srcw(w,Dt,D,src_est);
    w  = newton(fw,w,1e-6);
end
W  = oppKron2Lo(opKron(opDiag(w),opDirac(nsrc)),opDirac(nrec));

% evaluate misfit
[f,df]  = misfit(W*Dt - D);

% gradient
if nargout > 1
    g  = Jt'*(CC'*(W'*df));
end

% sanity check
if isdistributed(f)
    f = gather(f);
end

function [f,g,h] = srcw(w,Dt,D,fh)
% function handle for misfit as function of source weights

n2 = length(w);
n1 = length(Dt)/n2;

[f,g,h] = fh(oppKron2Lo(spdiags(w,0,n2,n2),opDirac(n1))*Dt-D);

g = oppKron2Lo(opDirac(n2),opOnes(1,n1))*(g.*conj(Dt));
h = oppKron2Lo(opDirac(n2),opOnes(1,n1))*(h.*conj(Dt).*Dt);

g = gather(g);
h = gather(h);


function [x] = newton(fh, x0, tol)
% Newton's method for problems with a strictly diagonal Hessian.
% 
% function [x] = newton(func, x0, tol)
%
% Input:  func  - function handle returning [value, gradient, hessian]
%                 where hessian is a vector containing the diagonal of the
%                 Hessian.
%         x0    - initial value of parameter
%         tol   - tolerance

% initialization
maxIter = 10;
mult    = 0.8;
c       = 0.01;
x       = x0;
iter    = 0;

[f, g, h] = fh(x);

while(norm(g) > tol)&&(iter < maxIter)
    
    % search direction
    d = -g./h;
    
    % Armijo linesearch
    step   = 1;
    lsiter = 0;
    while( fh(x + step*d) > f - c*step*abs(g'*d) )
        step = step*mult;
        lsiter = lsiter + 1;
        if lsiter > 20
            return;
        end
    end
    
    % update
    x = x + step*d;
    iter = iter + 1;
    [f, g, h] = fh(x);
 
end  % while
