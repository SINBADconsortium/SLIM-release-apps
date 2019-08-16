function x = EPSI_SLIM_pickSpikyFirstEvent(data, nt, nr, ns, OP_MULTPRED, OP_TIMEWIND, wavelet_length, options)        
% Pick the largest peak in every trace, so we get very sparse events for the first wavelet matching

        % Construct the operator invert in the initial spike-picking for first reflection event
        A = OP_MULTPRED * OP_TIMEWIND;
        
        % Pre-allocate a distributed empty solution if parallel processing is on
        if options.parallel
            preallocatedEmptySolutionVector = distVectorize(distributed.zeros(nr,ns,nt));
        else
            preallocatedEmptySolutionVector = [];
        end
        
        % Use a few iterations of SPGL1 to help condition the inital estimate (removes a lot of anti-causal events)
        % This obtains a rough sparse estimate of min |x|_1 s.t. Ax = data, where A is the multiple-prediction operator
        % 
        opts  = spgSetParms('iterations',        options.maxInitialIter, ...
                            'verbosity' ,        options.verbosity*2, ...
                            'bpTol'     ,     1e-5, ...
                            'decTol'    ,     5e-2, ...
                            'subspaceMin',    0,    ...
                            'quitPareto',     1,    ...
                            'ignorePErr',     options.spgl1_ignorePErr, ...
                            'nPrevVals' ,     4);
        x = spgl1(A, vectorizeData(data), [], [], preallocatedEmptySolutionVector, opts);
        
        % reshape solution into shot gathers
        if options.parallel
            x = distVec2distArray(x, [nr ns nt]);
            x = distSort_timeSlice2shot(x);
        else
            x = reshape(x, nt, nr, ns);
        end
        
        % assuming an impulsive spike event for the first reflection, pick the largest element per trace
        if options.parallel
            spmd
            for currShot = drange(1:ns)
                x(:,:,currShot) = pickLargestForEveryTrace(x(:,:,currShot));
            end
            end
        else
            for currShot = 1:ns
                x(:,:,currShot) = pickLargestForEveryTrace(x(:,:,currShot));
            end
        end
        
        % vectorize/resort after picking
        if options.parallel
            x = distSort_shot2timeSlice(x);
        end
        x = vectorizeData(x);        
        
        % Use exact line-search to find correct scaling for the picked spikes
        % Here we want to minimize the value of | b - sAx |_2^2, where b is the
        % recorded data, x is the picked first-event spikes, A is the multiple-
        % prediction operator OP_MULTPRED (remember here the Q term = 0), and s is
        % a scaling multiplier. The exact line-search assumes that the first-order
        % variational is 0, thus (Ax)'b = s(Ax)'Ax.
        
        A = OP_MULTPRED;
        b = vectorizeData(data);
        
        g = A' * b;
        y = A * x;
        s = (x' * g) / (y' * y);
        s_scaled = options.initialSpikeAmp * s;
        
        disp(['Initial picked spikes scaled by: ' num2str(undist(s)) ' * (' num2str(options.initialSpikeAmp) ' from options)'])
        disp(' ')
        
        x = s_scaled .* x;
        
end

function x = pickLargestForEveryTrace(x);
    % x should be a 2-d array, will pick out the number with largest abs value for every column
    % and zero the rest
    for currTrace = 1:size(x,2)
        % keep only the max amplitude value
        [maxVal, ind] = max(abs(x(:,currTrace)));
        maxVal_sign = sign(x(ind,currTrace));
        x(:,currTrace) = 0;
        x(ind,currTrace) = maxVal_sign * maxVal;
    end
end
     