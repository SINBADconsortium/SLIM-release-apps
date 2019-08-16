function [ L, R , res, time ] = spgLR( data, e, rank, sigma, opts )
%SPGLR - Distributed LR-matrix completion in the SPGL1 framework
% This function solves the following problem
%
%  min_{X}       nuclearNorm(X) 
%  subject to      ||A.*X - b|| <= sigma
%
% where
%   - A is a subsampling operator, represented as a sparse 0-1 matrix, specifying the location of the data
%   - b is the subsampled data, represented as a sparse matrix
%   - sigma is the noise level
%   - X is defined implicitly in terms of L and R factors, each with rank k, so that X = L * R'
% 
% using the SPGLR machinery detailed in 'FAST METHODS FOR DENOISING MATRIX COMPLETION FORMULATIONS, WITH APPLICATIONS TO ROBUST SEISMIC DATA INTERPOLATION', Aravkin, Kumar, et al. 
%
% Here L, R, b, and A are distributed arrays, so this function can be applied to very large matrices on distributed systems.
% 
% Usage: 
%    [ L, R ] = spgLR( data, e, rank, sigma, opts );
% 
% Input:
%    data  - nrows x ncols sparse matrix distributed along rows,
%            entries of the measured matrix
%    e     - 0/1 logical matrix, distributed along rows, 1 indicates
%            data at that point
%    rank  - rank of the factors
%    sigma - noise level to fit the data
%
%    opts  - struct with one or more of the following options
%          .maxIter    - maximum number of iterations (default: 500)
%          .Linit      - initial guess for the L factor (default: scaled gaussian noise)
%          .Rinit      - initial guess for the R factor (default: scaled gaussian noise)
%          .scale      - scale factor for random initialization, so each 
%                        entry is approximately this value (default: 1)
%          .verbosity  - silent (0), iteration(1, default), verbose(2)
%          .optTol     - optimality tolerance (default: 1e-4)
%          .decTol     - subproblem tolerance (default: 1e-6*opts.scale)
%          .M          - non-monotonicity parameter (default: 3)
%          .logFile    - log file to write output to (default: [])
%          .rankInc    - if lasso subproblem stalls away from boundary, the
%                        rank is increased by this amount (default: rank)
%          .maxStalls  - maximum number of times to increase rank when LASSO
%                        stalls, before giving up (default: 3)
%          .outputFreq - frequency with which to update output (default:
%                        every 50 iterations, or when a subproblem
%                        is solved)
%
% Output:
%   L, R - estimated LR factors, distributed
%
% Curt Da Silva, 2014
% curtd@math.ubc.ca
% 
    
assert(size(data,1) > 1 & size(data,2) > 1 & length(size(data))==2, 'Input data must be a matrix');
assert(size(e,1) == size(data,1) & size(e,2) == size(data,2) & length(size(e))==2 , 'Subsampling indices must be a logical matrix matching the size of the data matrix');

tic;    
LASSO_CONVERGED = 0; 
LASSO_UPDATE_TAU = 1;
LASSO_STALLED = 2;
EXIT_MAX_ITERATIONS = 3;

SILENT = 0; ITER = 1; VERBOSE = 2;

if exist('opts','var')==0 || isempty(opts), opts = struct; end
if isfield(opts,'maxIter'), maxIter = opts.maxIter; else maxIter = 500; end
if isfield(opts,'Linit'), L = opts.Linit; else L = []; end
if isfield(opts,'Rinit'), R = opts.Rinit; else R = []; end
if isfield(opts,'scale'), scaleFactor = opts.scale; else scaleFactor = 1; end
if isfield(opts,'verbosity'), verbosity = opts.verbosity; else verbosity = ITER; end    
if isfield(opts,'optTol'), optTol = opts.optTol; else optTol = 1e-4; end
if isfield(opts,'decTol'), decTol = opts.decTol; else decTol = 1e-6*scaleFactor; end
if isfield(opts,'M'), M = opts.M; else M = 3; end
if isfield(opts,'logFile'), logFile = opts.logFile; else logFile = []; end
if isfield(opts,'rankInc'), rankInc = opts.rankInc; else rankInc = rank; end
if isfield(opts,'maxStalls'), maxStalls = opts.maxStalls; else maxStalls = 3; end
if isfield(opts,'outputFreq'), outputFreq = opts.outputFreq; else outputFreq = 50; end
if isfield(opts,'penalty'), penalty = opts.penalty; else penalty = @LSPenalty; end
if isfield(opts,'stagCheck'),stagCheck = opts.stagCheck; else stagCheck = 10; end

nlabs = parpool_size(); 
assert( nlabs > 0, 'This function can only be run with an open parallel pool');

nrows = size(data,1); ncols = size(data,2);
rowidx = ceil(linspace(1,nrows+1,nlabs+1));
colidx = ceil(linspace(1,ncols+1,nlabs+1));
stepMin = 1e-10;
stepMax = 1e4;
theta = 0.25;

if ~isempty(logFile)
    fid = fopen(logFile,'w');
else
    fid = -1;
end
nStalls = 0;
fevals = 0;
gevals = 0;
bNorm = norm(data,'fro'); 

output(['------------------------ SPGLR started on ' datestr(now) ' ------------------ '],ITER);
output(['Maximum iterations ' num2str(maxIter)],ITER);
output(['Data norm: ' num2str(bNorm,'%3.3e')],VERBOSE);
output(['Optimality tolerance: ' num2str(optTol,'%3.3e')],VERBOSE);
output(['Subproblem tolerance: ' num2str(decTol,'%3.3e')],ITER);
output(['Minimum step: ' num2str(stepMin,'%3.3e')],VERBOSE);
output(['Maximum step: ' num2str(stepMax,'%3.3e')],VERBOSE);

% Initialization of L, R factors
spmd
    minRow = rowidx(labindex); maxRow = rowidx(labindex+1)-1;
    minCol = colidx(labindex); maxCol = colidx(labindex+1)-1;
    Lloc = sqrt(scaleFactor)*randn(maxRow-minRow+1,rank)/sqrt(rank);         
    L2 = codistributed.build(Lloc,codistributor1d(1,diff(rowidx),[nrows,rank]),'noCommunication');  
    
    if isempty(L)
        L = L2;
    else
        L = redistribute(L,getCodistributor(L2));
        L2 = [];
    end
    Rloc = sqrt(scaleFactor)*randn(maxCol-minCol+1,rank)/sqrt(rank); 
    R2 = codistributed.build(Rloc,codistributor1d(1,diff(colidx),[ncols,rank]),'noCommunication');  
    if isempty(R)
        R = R2;
    else
        R = redistribute(R,getCodistributor(R2));
        R2 = [];
    end    
    % Ensure input has the correct distribution
    Lcodist = getCodistributor(L);
    
    data = redistribute(data,codistributor1d(1,Lcodist.Partition,[nrows,ncols]));
    dataLoc = getLocalPart(data); 
    e = redistribute(e,codistributor1d(1,Lcodist.Partition,[nrows,ncols]));
    eloc = getLocalPart(e);
    nnzloc = length(find(eloc));    
end

%tau = 0.5 * (nrows + ncols) * scaleFactor;
%[L,R] = project(L,R,tau);
tau = 0.5*(norm(L,'fro')^2+norm(R,'fro')^2);
nnz = zeros(nlabs,1); 
for j=1:nlabs
    nnz(j) = nnzloc{j}; 
end       

iter = 0;

    function output(msg,level)
        if verbosity >= level
            disp(msg);
            if fid > 0
                fprintf(fid,[msg '\n']);
            end
        end
    end
     
    function t = innerprod(L1,L2)
        spmd,
            loc1 = getLocalPart(L1); loc2 = getLocalPart(L2); tloc = vec(loc1)' * vec(loc2);
            t = gplus(tloc,1);
        end
        t = t{1};        
    end
   
    
    function [L,R] = project(L,R,tau)
        n = sqrt(norm(L,'fro')^2+norm(R,'fro')^2);
        if n > sqrt(2*tau)
           L = L*sqrt(2*tau)/n; R = R*sqrt(2*tau)/n;
        end
    end    
    
    function [L,R,f,lassoErr,locIter] = lassoSubproblem(L,R,tau)
        locIter = 0;
        testUpdateTau = false;
        
        %Initial function evaluation
        [f,gL,gR] = spgLRobj(L,R,data, e, penalty);       
        fevals = fevals+1; gevals = gevals+1;
        gNorm = sqrt(norm(gL,'fro')^2+norm(gR,'fro')^2);
        if gNorm < 1/stepMax
            gStep = stepMax;
        else
            gStep = min( stepMax, max(stepMin,1/gNorm) );
        end
        %history of function values
        lastFv = -inf(M,1);
        lastFv(1) = f; 
        fOld = f; fMin = f;
        fStag = f;
        
        lassoErr = -inf;
        optGap = [];
        output(sprintf('itr % 5d f: %6.5e g:    %5.2e  tau: %3.3e',iter+locIter,f,sqrt(norm(gL,'fro')^2+norm(gR,'fro')^2),tau),ITER);   
        while 1
            locIter = locIter + 1;
            rNorm = f;            
            rErr = (rNorm-sigma)/max(1,rNorm);
            
            testRelChange1 = (abs(f-fOld) <= 1e-1*decTol * f);
            testRelChange2 = (abs(f-fOld) <= 1e-1 * f * (abs(rNorm-sigma)));
            testUpdateTau = ((testRelChange1 && rNorm > 2*sigma) || ...
                (testRelChange2 && rNorm <= 2*sigma)) && ~testUpdateTau && locIter > 1;
            testDone = rErr <= optTol || rNorm < optTol * bNorm;
            
            % Stagnation check
            if locIter > 1 && mod(locIter,stagCheck)==0
                if abs(f-fStag) <= 1e-2 * abs(fStag)
                    lassoErr = LASSO_UPDATE_TAU; break; 
                else
                    fStag = f;
                end
            end
            
            if testDone 
               lassoErr = LASSO_CONVERGED;
               output('Lasso converged',ITER);
               break; 
            end            
            
            if testUpdateTau
               Lnew = L-gL; Rnew = R-gR;
               [Lnew,Rnew] = project(Lnew,Rnew,tau);
               Lnew = Lnew-L; Rnew = Rnew-R;
               optGap = gather(max( max(max(abs(Lnew))), max(max(abs(Rnew)))));
               if optGap < decTol
                   if 0.5*(norm(L,'fro')^2+norm(R,'fro')^2) < 0.99*tau
                       lassoErr = LASSO_CONVERGED;
                       output('Optimality gap reached at interior of feasible set, further progress not achieved, quitting',ITER);
                       break;
                   end
                   lassoErr = LASSO_UPDATE_TAU; break; 
               end              
            else
                optGap = [];
            end
            
            % Non monotone line search
            fMax = max(lastFv);            
            lsIter = 1; maxLSiter = ceil((log(stepMax)-log(stepMin))/log(1/theta)); lserr = 0; 
                        
            while 1
                if lsIter > maxLSiter
                    lserr = inf;
                    break; 
                end
                Lnew = L - gStep * gL; Rnew = R - gStep * gR; 
                n = sqrt(norm(Lnew,'fro')^2+norm(Rnew,'fro')^2);
                if n > sqrt(2*tau)
                    Lnew = Lnew*sqrt(2*tau)/n; Rnew = Rnew*sqrt(2*tau)/n;
                end
                s = Lnew-L; gts = innerprod(s,gL); s = Rnew-R; gts = (gts + innerprod(s,gR));
                
                fNew = spgLRobj(Lnew,Rnew,data, e, penalty);
                fevals = fevals+1;
                gamma = 1e-4;
                
                if n > sqrt(2*tau)                    
                    output(sprintf('LSitr % 3d f: %6.5e step: %1.2e  tau: %2.3e',lsIter,fNew,gStep,tau),VERBOSE)
                else
                    output(sprintf('LSitr % 3d f: %6.5e step: %1.2e  tau: %2.3e',lsIter,fNew,gStep,0.5*n^2),VERBOSE);                        
                end
                
                if fNew-fMax < gamma*gStep*gts
                    break; 
                end
                gStep = gStep * theta;
                lsIter = lsIter + 1;
            end                            
            
            if lserr
                %Nothing fancy at the moment
                disp(['Line search failed with error ' num2str(lserr) ]);                
                break;
            else
                %No line search issues, update L + R
                fOld = f; 
                gLold = gL; gRold = gR;                
                
                [f,gL,gR] = spgLRobj(Lnew,Rnew,data, e, penalty);
                fevals = fevals+1; gevals = gevals+1;
                lastFv(mod(locIter,M)+1) = f;
                fMin = min(f,fMin);     
                if(mod(iter+locIter,outputFreq)==0)
                    if isempty(optGap)                        
                        output(sprintf('itr % 5d f: %6.5e ',iter+locIter,f),ITER);                           
                    else
                        output(sprintf('itr % 5d f: %6.5e opt:  %3.2e',iter+locIter,f,optGap),ITER);           
                    end
                end
                s = Lnew - L; sts = innerprod(s,s);
                y = gL - gLold; sty = innerprod(s,y);  
              
                s = Rnew - R; y = gR - gRold; sts = sts + innerprod(s,s); sty = sty + innerprod(s,y);  
              
                clear s y Lold gLold;

                L = Lnew; R = Rnew;
            end                       
                       
            if sty <= 0, 
                gStep = stepMax; 
                output('Neg curvature',VERBOSE);
            else gStep = min( stepMax, max(stepMin, sts / sty ) );
            end           
            
            if iter+locIter >= maxIter
               lassoErr = EXIT_MAX_ITERATIONS;
               break;
            end
        end        
     
    end    
    
    output(['Target data misfit ' num2str(sigma,'%3.3e')],ITER);
    output(['Initial tau ' num2str(tau,4)],ITER);
    
    while 1                        
        
        [L,R,fk,err,locIter] = lassoSubproblem(L,R,tau);
        iter = iter + locIter;
        if err == LASSO_UPDATE_TAU
            nStalls = 0;
            tauOld = tau;  
            output('Tau update',ITER); 
            h = 1e-8; dtau = max(h*tau,h);
            [L,R,fnew,~,locIter] = lassoSubproblem(L,R,tau+dtau);

            tau = max(0, tau - (fk-sigma)/((fnew-fk)/(dtau))); %Secant approximation to update tau
            output(['New tau : ' num2str(tau,4)],ITER);
            iter = iter + locIter;            
        elseif err == LASSO_STALLED
            if nStalls >= maxStalls
                output('Rank increases stalled too many times, quitting',SILENT);
                break;
            else
                % Before stalling too much, try to increase the
                % size of L and R. EXPERIMENTAL
                output('Lasso stalled',ITER);
                nStalls = nStalls + 1;
                
                % Assert that the L, R factors are full rank
                irows = sort(randperm(nrows,round(2*rank*log(rank))),'ascend');                
                Lsub = gather(L(irows,:));
                sv = svd(Lsub); effrank = length(find(sv >= 1e-6*sv(1)));
                
                % If they're not, then exit because something has gone horribly wrong
                if effrank < rank
                    output('Increasing rank wont decrease objective, quitting ',SILENT);
                    break;
                end
                output(['Increasing rank after ' num2str(nStalls) ' stalls'],ITER);                                                                                                
                spmd,
                    Lloc = getLocalPart(L); 
                    Lloc = [Lloc,1e-3 * sqrt(scaleFactor)* randn(size(Lloc,1),rankInc)/sqrt(rankInc)];
                    L = codistributed.build(Lloc,codistributor1d(1,diff(rowidx),[nrows,rank+rankInc]),'noCommunication');
                    Rloc = getLocalPart(R);
                    Rloc = [Rloc,1e-3 * sqrt(scaleFactor)* randn(size(Rloc,1),rankInc)/sqrt(rankInc)];
                    R = codistributed.build(Rloc,codistributor1d(1,diff(colidx),[ncols,rank+rankInc]),'noCommunication');
                end
                
                rank = rank + rankInc;
                [L,R] = project(L,R,tau);
            end
        elseif err == LASSO_CONVERGED
            output(['Converged'],ITER);
            break; 
        elseif err == EXIT_MAX_ITERATIONS
            output('Too many iterations',ITER);
            break;
        else
           output(['Error ' num2str(err) ' from lasso subproblem'],SILENT); 
           break;
        end
    end
    res = fk;
    time = toc;
    output(['Total time ' num2str(time/3600) ' hours for ' num2str(iter) ' iterations'],ITER);
    if fid > 0
        fclose(fid);
    end
    
end

