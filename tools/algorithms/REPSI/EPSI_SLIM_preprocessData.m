function [data, nt, nr, ns, options] = EPSI_SLIM_preprocess_Data(data, dt, options)
    
    % for parallel processing, make 'data' a native distributed array
    if options.parallel
        disp('The use of Matlab Parallel Computing Toolbox is ON')
        
        % test to see if MATLAB PCT exists and workers are active
        try
            poolSize = parpool_size();
            if poolSize == 0
               error('RobustEPSI:parallelExecution:noWorkers','There are no parallel pool workers abailable, please make sure a valid parallel pool session is active before attempting parallel execution') 
            else
                disp(['... found ' int2str(poolSize) ' parallel pool workers'])
                disp(' ')
                % disable unnecessary warnings
                warning off distcomp:spmd:RemoteTransfer 
            end
        catch ME
            error('RobustEPSI:parallelExecution:noPCT','Command MATLABPOOL does not exist, this installation of MATLAB does not have access to Parallel Computation Toolbox')
        end
        
        data = distributed(data);
    end
    
    % get data dimensions
    dims = size(data);
    nt   = dims(1);
    nr   = dims(2);
    ns   = dims(3);
    
    % make sure the data is in correct precision, also set precision flag
    if options.singlePrecision
        data = single(data);
        options.precision = 'single';
        % projection step in spgl1 might be imprecise, choose to allow this
        options.spgl1_ignorePErr = 1;
    else
        data = double(data);
        options.precision = 'double';
        options.spgl1_ignorePErr = 1;
    end

    % renormalize the data if desired
    if options.renormData
        datamax = norm(vectorizeData(data),inf);
        data = data ./ datamax;
    end

    % due to Wavelet sparsity transform needing powers of 2, we sometimes need to pad nt to the next largest power of 2
    if options.padtime ~= 0
        if (options.padtime > 0)  % pad
            data = [data; zeros(options.padtime,nr,ns,options.precision)];
        else
            data = data(1:(nt+options.padtime),:,:); % chop
        end
        nt = nt+options.padtime;
    end

    % Shift data downwards if needed (can fake anti-causal part of source this way)
    if options.downshift > 0
        data = [zeros(options.downshift,nr,ns,options.precision); data];
        data = data(1:end-options.downshift,:,:);
    end
    
    % Weight the traces as a power of time
    if options.timeweightgamma
        data = timeweighting(data,options.timeweightgamma,dt);
    end

    % Kill the near-offsets if reqired
    if (options.num_traces_missing > 0)
        for num_shotindex = 1:ns
            list_tracepos = [(num_shotindex - options.num_traces_missing) : (num_shotindex + options.num_traces_missing)];
            data(:,list_tracepos(list_tracepos>0 & list_tracepos<=ns),num_shotindex) = 0;
        end
    end
    
    % For parallel processing, re-sort data into time-slice gathers for the rest of the processing
    if options.parallel
        data = distSort_shot2timeSlice(data);
    end

end
