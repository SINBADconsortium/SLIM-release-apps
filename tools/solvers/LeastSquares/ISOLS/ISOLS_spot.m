function [Uk,iter,resvec,timing,residual,W,Y]=ISOLS_spot(Q,D,lambda,P,H,params)
% ISOLS solves least-squares problems of the form:
% argmin_u || [H ; lambda*P]u - [Q;lambda*D] ||_2
% using iterative methods implemented though the SPOT linear operator
% toolbox.
%
% Author: Bas Peters
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: August 2015.
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
% If you have any questions, errors or disappointing results, email
% (bpeters {at} eos.ubc.ca)
%
% Input:
%           Q       - source matrix. size(Q,1) must match source grid
%                     definition, size(Q,2) determines the number of
%                     sources.
%           D       - Observed data matrix (nsrc X nrec)
%           lambda  - penalty parameter (scalar, >0 )
%           P       - Receiver (observation) matrix, (nrec X N_grid_points)
%           H       - SPOT operator, action of which defines H^(-1) and H^(-*) on a
%                     vector using an iterative method.
%           params  - structure with some options
%           params.analysis  - output residual vectors or not
%           params.sim_rec   - use simultaneous receivers
%           params.n_sim_rec - number of simultaneous receivers to use
%           params.refine    - use solution refinement at the end
%           params.maxit_refine - number of refinement iterations

% Notes: 
% -currently only outputs the field Uk and timing other parameters
%  are dummy (will be fixed later)
% -ISOLS_spot : iterative method for slightly overdetermined systems, SPOT
% linear operator toolbox version.

% Parse settings
resvec    = 0;
sim_rec   = 0;
n_sim_rec = 0;
if exist('params','var') && ~isempty(params)
    if isfield(params,'sim_rec'),  sim_rec  = params.sim_rec; end
    if isfield(params,'n_sim_rec'),  n_sim_rec  = params.n_sim_rec; end
end

%% initialize
residual= 0;
iter.u  = 0;
iter.w  = 0;
tic;    %timer for 'common part' per frequency

%% setup up randomization and subsampling

%currently only random subsampling is implemented (in the form of a random
%selection of the receiver/data blocks after multiplying by DCT*random[+1
%-1]
if sim_rec==1
    nrec = size(P,1);
    list_integers = [ones(ceil(nrec/2),1); -1*ones(floor(nrec/2),1)];
    ind = randperm(length(list_integers),length(list_integers));
    Dr = diag(list_integers(ind));
    F = opDCT(nrec);
    R = eye(nrec);
    ind=randperm(nrec,n_sim_rec);
    R = R(ind,:);
    V = R*F*Dr;
    
    P=V*P;
    D=V*D;
end



%% compute H^-*P explicitly
for i=1:size(P,1)
    W(:,i) =H'\(full(P(i,:)'));
end
iter.w=0;%iter.w+length(resvec.W);

%% Sherman-Morrison-Woodbury part of the algorithm
X = W'*W;
S = inv(speye(size(W,2))+ lambda^(-2)*X);
Y = lambda*Q+lambda^(-1)*W*(S*(-W'*Q+D));

timing.common=toc;

%% for all right hand sides
tic; %timer for right hand sides part

%obtain fields from y
for i=1:size(Q,2)
    i
    Uk(:,i)= H\Y(:,i);
end
iter.u=0;
Uk=(1/lambda).*Uk;
timing.rhs=toc;


