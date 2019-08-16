function [dm,info] = linbregLSRTM(dData,model,m0,q,precon,options)
%linbregLSRTM Sparsity promoting LSRTM in the time domain using the
% linearized Bregman method.
%
%   [dm,info] = linbregLSRTM(data,model,m0,q,precon,options)
%
%--------------------------------------------------------------------------
%
% INPUT:
%   data     - input shot record(s) of dimensions ns x nrec x nsrc (ideally 
%              this is linearized, i.e. single scattered data)
%
%   model    - structure with model parameters
%
%   m0       - smooth background model (squared slowness [s^2/km^2])
%
%   q        - source wavelet
%
%   precon   - structure of preconditioners (set to '1' to enable)
%              .modelTopmute      model space topmute
%              .modelDepth        model space depth scaling
%              .dataScaling       date space half integration
%              .dataTopmute       data space topmute
%
%   options  - structure for LSRTM options
%              .iterations        number of iterations (default is 10)
%              .lambda            soft thresholding parameter (default)
%              .sigma             L2 norm of noise (default is 0)
%              .numShots          shots per iteration (default is 1)
%              .snapshotX         path to snapshots of primal variable
%              .snapshotZ         path to snapshots of dual variable
%              .dm                true image to calculate model error
%
%
% OUTPUT:
%   dm       - model perturbation/image
%   info     - structure with information (no. of iterations, lambda,
%              relative residual, model error)
%
%
%
% Author: Philipp Witte
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: February, 2016

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.


% Data dimensions
nsrc = length(model.xsrc);
nrec = length(model.xrec);
ns = length(model.NyqT);

% Check input options and preconditioners and set default values
[options, precon] = checkOptions(options,precon);

% Snapshot operators
if isfield(options,'snapshotZ')
    opSaveZ = opSaveSnapshot(model.n(1)*model.n(2),options.snapshotZ);
end

if isfield(options,'snapshotX')
    opSaveX = opSaveSnapshot(model.n(1)*model.n(2),options.snapshotX);
end

%% Preconditioners for LSRTM

N = prod(model.n);

% Depth preconditioner
if precon.modelDepth 
    Depth = opDiag(sqrt(model.d(1):model.d(1):model.d(1)*model.n(1)));
    opDepth = opKron(opDirac(model.n(2)),Depth);
else
    opDepth = opDirac(N);
end

% Model topmute
if precon.modelTopmute
    idxWaterCol = model.zrec(1)/model.d(1);
    taper = 20;    % no. of points for taper
    mute = opLinMute(model.n(1),idxWaterCol-taper,idxWaterCol);
    opMute = opKron(opDirac(model.n(2)),mute);
else
    opMute = opDirac(N);
end

% Right hand preconditioner
opRight = opMute*opDepth;

% Left hand preconditioner (data space), half integration
if precon.dataScaling
    H = opHalfInt(model);
else
    H = opDirac(ns*nrec);
end

% Set up curvelet transform
NBscales = max(1,ceil(log2(min(model.n(1:2))) - 3));
NBangles = 16;
Finest = 1;
Ttype = 'ME';
IS_real = 1;
C = opCurvelet(model.n(1),model.n(2),NBscales,NBangles,Finest,Ttype,IS_real);


%% Linearized Bregman

% Some initializations 
x = zeros(length(m0(:)),1);
z = x;
numTotalShots = nsrc;
numSubShots = options.numShots;
modelSub=model;
sigma = options.sigma*numSubShots/numTotalShots;
H = opKron(opDirac(numSubShots),H);

% Shuffle the shot records
src = randperm(numTotalShots);

% soft thresholding function
softThreshold = @(x,lambda) max(abs(x)-lambda,0).*sign(x);

disp('Start LSRTM')
for i=1:options.iterations
    disp(i)
    
    % Subsample shots
    idx0 = mod(i,numTotalShots);
    if idx0==0
        idx0 = numTotalShots;
    end
    idx1 = idx0*numSubShots-numSubShots+1;
    idx2 = idx1+numSubShots-1;
    modelSub.xsrc = model.xsrc(src(idx1:idx2));
    modelSub.ysrc = model.ysrc(src(idx1:idx2));
    modelSub.zsrc = model.zsrc(src(idx1:idx2));
    
    % Subsample born operator and top mute operator
    J = opBorn(m0,modelSub,q);
    if precon.dataTopmute
        Td = opTopMute(modelSub);
    else
        Td = opDirac(ns*nrec*numSubShots);
    end
    
    % Composite matrix with all preconditioners
    opLeft = Td*H;
    J = opLeft*J*opRight;
    
    % Subsample input data and apply left-hand preconditioning
    b =opLeft*vec(dData(:,:,src(idx1:idx2)));
  
    % Demigration J*x to calculate residual r
    if i==1
        % skip first demigration since initial guess is zero
        r = -b;  
    else
        r = J*x-b;
        % Alternatively: r = opLeft*forwardModel(m0+opRight*x,model,q)
    end
    
    % Migration of residual: J'*(J*x-b)
    JTr = J'*r;
    % Alternatively: JTr = opRight'*RTM(m0,model,q,opLeft'*r)

    % Stepsize
    t = norm(r,2)^2/norm(JTr,2)^2;
    
    % determine lambda if not supplied
    if i==1 && isempty(options.lambda)
        lambda = 0.1*max(abs(C*(t*JTr*max(0,1-sigma/norm(r,2)))));
    end
    
    % Update variables
    z = z - t*JTr*max(0,1-sigma/norm(r,2));
    x = C'*softThreshold((C*z),lambda);

    % Save snapshots
    if isfield(options,'snapshotX')
        opSaveX*opRight*x;
    end
    if isfield(options,'snapshotZ')
        opSaveZ*opRight*z;
    end
 
    % Current residual and error
    info.residual(i) = norm(r,2);
    if isfield(options,'dm')
        info.err(i) = norm(x-options.dm,2);
    end
    
end


%% Gather info and reshape result
info.iterations = options.iterations;
info.lambda = lambda;
info.sources = src;
dm = reshape(opRight*x,model.n); 


end



function [options,precon] = checkOptions(options,precon)
% CHECKOPTIONS Set default options if not supplied
%

% 10 iterations as default
if ~isfield(options,'iterations')
    options.iterations = 10;
end

% automatic lambda determination in first iteration as default
if ~isfield(options,'lambda')
    options.lambda = [];
end

% zero noise as default
if ~isfield(options,'sigma')
    options.sigma = 0;
end

% one shot per iteration as default
if ~isfield(options,'numShots')
    options.numShots = 1;
end

% disable all preconditioners as default
if ~isfield(precon,'modelTopmute')
    precon.modelTopmute = 0;
end

if ~isfield(precon,'dataTopmute')
    precon.dataTopmute = 0;
end

if ~isfield(precon,'modelDepth')
    precon.modelDepth = 0;
end

if ~isfield(precon,'dataScaling')
    precon.dataScaling = 0;
end


end


