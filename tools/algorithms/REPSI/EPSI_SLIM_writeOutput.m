function EPSI_SLIM_writeOutput( opEPSI, data, xphys, q_est, q_total, nt, nr, ns, dt, wavelet_window_start, wavelet_window_end, total_itercount, residual_log, options )
% Prepare the output of EPSI_SLIM and write to file
    
    disp(' ')
    
    % prepare Green's function
    if options.parallel
        primaryIR = distVec2distArray(xphys, [nr ns nt]);
    else    
        primaryIR = reshape(xphys, [nt,nr,ns]);
    end
    clear xphys
    
    if (options.timeweightgamma ~= 0)
        primaryIR = timeweighting(primaryIR,options.timeweightgamma,dt,'inv');
    end
    
    [OP_EPSI OP_PTERM OP_QTERM OP_XTERM] = opEPSI(primaryIR,q_est(:,end),0,dt,options.useOblique);
    
    % failsafe
    % if options.savesolmat
    %     save(options.sol_file, 'primaryIR', 'q_est', 'q_total');
    % end
    % if options.savepreviewmat
    %     preview_primaryIR = primaryIR(:,:,ceil(ns/2));
    %     save(options.preview_file, 'preview_primaryIR', 'q_est', 'q_total');
    % end
    
    % Define how to render final output
    primary = OP_QTERM * vectorizeData(primaryIR);
    
    % distributed version requires a reshape and re-sort/permute
    if options.parallel
        primary = distVec2distArray(primary, [nr, ns, nt]);
        
        primaryIR = distSort_timeSlice2shot(primaryIR);
        primary = distSort_timeSlice2shot(primary);
    else
        primary = reshape(primary,nt,nr,ns);
    end
    
    % change variable name for recorded data, and prepare for output
    initial_data = data;
    clear data

    if (options.timeweightgamma ~= 0)
        initial_data = timeweighting(initial_data,options.timeweightgamma,dt,'inv');
    end
    
    if options.parallel
        initial_data = distSort_timeSlice2shot(initial_data);
    end
    
    %% For debugging output, save all the relavant solution info to matlab file
    if options.savepreviewmat
        preview_primary = undist(primary(:,:,ceil(ns/2)));
        preview_primaryIR = undist(primaryIR(:,:,ceil(ns/2)));
        preview_initialdata = undist(initial_data(:,:,ceil(ns/2)));
        zoffset_primary = zeros(nt,ns,options.precision);
        zoffset_initialdata = zeros(nt,ns,options.precision);
        for k = 1:ns
            zoffset_primary(:,k) = undist(primary(:,k,k));
            zoffset_initialdata(:,k) = undist(initial_data(:,k,k));
        end
        
        save(options.preview_file,'zoffset_initialdata','zoffset_primary','preview_primary','preview_primaryIR','preview_initialdata','q_est','q_total','options','dt','nt','nr','ns','wavelet_window_start', 'wavelet_window_end','total_itercount','residual_log');
        disp('Recording solution information... preview mat file saved successfully')

        if options.show_preview
            disp('    (preview data plotted to screen)')
            display_Preview(options.preview_file)
        end
    end

    if options.savesolmat
        % warning, will gather distributed solution to the head node
        primary = undist(primary);
        primaryIR = undist(primaryIR);
        initial_data = undist(initial_data);
        save(options.sol_file,'primary','primaryIR','q_est','q_total','initial_data','options','dt','wavelet_window_start', 'wavelet_window_end','total_itercount','residual_log');
    	disp('Recording solution information... solution mat file saved successfully')
    end

    %% Write to output files:
    
    % final result (primary wavefield)
    if ~isempty(options.output_primary_file)
        primary = undist(primary);
        write_EPSIdata(options.output_primary_file, primary, dt)
        clear primary
    end
    
    % final result (primary impulse response wavefield, deconvolved from the source signature)
    if ~isempty(options.output_primaryIR_file)
        primaryIR = undist(primaryIR);
        write_EPSIdata(options.output_primaryIR_file, primaryIR, dt)
        clear primaryIR
    end
    
    % estimated wavelet
    if ~isempty(options.output_wavelet_file)
        wavelet = q_est(:,end);
        
        if options.useOblique
            % correct for obliquity factor
            nt_conv = length(wavelet);
            F = opFFTsym_conv(nt_conv);
            nf = size(F,1);
            Binv = opObliqInv(dt,nf,1);
            wavelet_invObliq = F' * Binv * F * wavelet;
                
            write_EPSIwavelet(options.output_wavelet_file, wavelet, nt, dt, wavelet_window_start, wavelet_window_end, wavelet_invObliq)
        else
            write_EPSIwavelet(options.output_wavelet_file, wavelet, nt, dt, wavelet_window_start, wavelet_window_end)
        end
    end

end