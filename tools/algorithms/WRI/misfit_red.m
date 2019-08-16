function [f,g,h,w,f_aux] = misfit_red(m,Q,D,model,params)

% misfit for FWI with source estimation
%
% use:
%   [f,g,w] = misfit_red(m,Q,D,model,params)
%
% input:
%   m  - reference model
%   dm - model perturbation
%   Q  - source
%   D  - data  
%   model - struct with model paramaters
%   params.C         - aqcuisition mask as matrix of size nrec x nsrc
%   params.srcest    - if true, use source estimation (default: false)
%   params.computeLU - use LU factors for inverting helmholtz (default: false)
%   params.nthreads  - # of threads to use in spmd blocks (default: 0, single threading)
% output:
%   f - function value
%   g - gradient
%   w - estimated source wavelet
%
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: 2012
%  
% Updated by: Curt Da Silva, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


nsrc  = size(Q,2);
nrec  = length(model.xrec)*length(model.zrec);
nfreq = length(model.freq);
if exist('params','var') && isfield(params,'C'), C = params.C; else, C = ones(nrec,nsrc); end
if exist('params','var') && isfield(params,'srcest'),srcest = params.srcest; else, srcest = false; end

% model data
[Dt,Jt] = F(m,Q,model,params);

D  = reshape(D, [nrec*nsrc,nfreq]);
Dt = reshape(Dt,[nrec*nsrc,nfreq]);

% source estimation
w  = ones(1,nfreq);
f  = zeros(1,nfreq);
dR = D;
for k = 1:nfreq
    at = gather(C(:).*Dt(:,k));
    a  = gather(C(:).*D(:,k));
    if srcest
        fhk  = @(w)srcest(w,at,a,@(x)twonorms(x));
        w(k) = newton(fhk,1,1e-6);
    end
    [f(k),dfk] = twonorms(w(k)*at - a);
    dR(:,k)    = conj(w(k))*C(:).*dfk;
end

f = sum(f);

if nargout > 1        
    g = Jt'*vec(dR);
end

%initialize some dummy variables
f_aux.dat=0;
f_aux.pde=0;
h     = oppH(m,Q,D,model,params);

end

function [f,g,h] = srcest(w,a,b,fh)
    [f,df,ddf] = fh(w*a - b);
    f = gather(f);
    g = gather(a'*df);
    h = gather(a'*(ddf.*a));
end

