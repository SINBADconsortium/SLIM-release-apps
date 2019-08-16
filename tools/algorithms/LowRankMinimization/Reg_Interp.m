function [ output] = Reg_Interp(A,para,options,rank,sigmafact,initfact,fid)
% This function perform regularization and interpolation on monochromatic frequency slices using 
% Generalized SPGL1.
%
% use:
%   [output] = Reg_Interp(A,para,opts,rank,sigmafact,initfact)
%
% input:
% A         - input data at irregular grid with or without missing shots/receivers.
% para      - This will contain the information about data.
% opts      - parameters for SPGL1
% rank      - rank of the system to be used in the processing. It can be a
%             single number or a vector equal to the number of frequency.
% sigmafact - defines the percentage of noise tolerance in the output
%             results
% initfact  - This factor control the amplitude of intial random guess. it 
%             is used to make the amplitude of intial guess to be equal 
%             to the amplitude of input data. 
% freqindx  - save the output log file from SPGL1 for each frequency
% output:
%   Data    - Regularized and Interpolated data

% Author: Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date : April, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if para.penalty==1
    params.nu = para.nu; % This value define the amplitude of noise present in the data.
end

nR                = para.nrow;
nC                = para.ncol;
MH                = opMH(nR,nC);
params.numr       = nR;
params.numc       = 2*nC-1;
D                 = vec(A);
params.Ind        = find(D==0);
params.afun       = @(x)afun(para.Fireg*MH'*vec(x),params.Ind); % A
params.afunT      = @(x)reshape(MH*para.Fireg'*vec(x),nR,2*nC-1); % A'
params.funForward = @NLfunForward;
% Run the main algorithm
params.nr = rank;
LInit     = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit     = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit     = initfact*[vec(LInit);vec(RInit)]; % Initial guess
tau       = norm(xinit,1);
sigma     = sigmafact*norm(D,2);
opts = spgSetParms('optTol',options.optTol, ...
                   'bpTol', options.bpTol,...
                   'decTol',options.decTol,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', options.iteration,...
                   'fid',fid,...
                   'weights', []);
if para.penalty==0
    opts.funPenalty = @funLS;
else
    opts.funPenalty = @funST;
end
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XR        = spgl1(@NLfunForward,D,tau,sigma,xinit,opts,params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e         = params.numr*params.nr;
L1        = XR(1:e);
R1        = XR(e+1:end);
L1        = reshape(L1,params.numr,params.nr);
R1        = reshape(R1,params.numc,params.nr);
output    = reshape(para.Freg*MH'*vec(L1*R1'),nR,nC);
end

