function [D,J] = F(m,Q,model,params,freqsxsy)
% Frequency domain FD modeling operator
%
% use: 
%   [D,J] = F(m,Q,model,{params})
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source and has to be
%                       distributed over the last dimension.
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%
%   params            - optional struct of performance options
%   params.computeLU  - if true, will use LU factors for inverting Helmholtz (default: false)
%   params.nthreads   - if > 0, forces matlab to use multithreading inside spmd blocks (default: 0)
%   params.helm_save_dir - if not empty and params.computeLU == true, then
%                          specifies the directory where the helmholtz
%                          factorizations will be stored/loaded,
%                          corresponding to the model m
% output:
%   D  - Data cube (nrec x nsrc x nfreq) as (distributed) vector. nsrc  = size(Q,2);
%                                                                 nrec  = length(zrec)*length(xrec) 
%                                                                 nfreq = length(freq)
%   J  - Jacobian as SPOT/pSPOT operator
%%
% *Example*
% Below example defines a grid [0,1000]^2 with 10 m gridspacing. The
% velocity is 2 km/s. The sources and receivers coincide at x = [0:10:1000]
% and z = 10m.
%
%%
% model.o = [0 0];
% model.d = [10 10];
% model.n = [101 101];        
% model.nb = [10 10];
% model.freq = [5 10 15 20]; 
% model.f0 = 0;
% model.zsrc = 10; 
% model.xsrc = 0:10:1000;
% model.zrec = 10; 
% model.xrec = 0:10:1000;
% m = .25*ones(prod(model.n),1);
% Q = speye(length(model.xsrc));
% D = F(m,Q,model);
% 
%
% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% Updated by: Curt Da Silva, Bas Peters, 2015
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if exist('params','var')==0, params = struct; end
params.wri = false;
params = default_fwi_params2d(params);
model = fwi2d_model_compatibility(model);

if exist('freqsxsy','var')==0, freqsxsy = []; end

nlabs = parpool_size();
if nlabs==0
    J = opDF(m,Q,model,params);
    D = PDEfunc(PDEopts.FORW_MODEL,m,Q,[],[],model,params,freqsxsy);
else    
    J = oppDF(m,Q,model,params);            
    D = PDEfunc_dist(PDEopts.FORW_MODEL,m,Q,[],[],model,params,freqsxsy);    
end
D = vec(D);

