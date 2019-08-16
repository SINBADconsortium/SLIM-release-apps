function display_Preview(filename)
% Reads an EPSI preview matlab file and plot some results

load(filename,'zoffset_initialdata','zoffset_primary','preview_primary','preview_primaryIR','preview_initialdata','q_est','q_total','options','dt','wavelet_window_start', 'wavelet_window_end');

nt = size(preview_primaryIR,1);
time_tick_interval = 0.5;
time_tick_persamp = floor(time_tick_interval/dt);
time_axis = [1:time_tick_persamp:nt];
time_axis_label = [0.0:0.5:nt*dt]';


%% Plot wavelet info

% Preprocess wavelet, window it out
wavelet = q_est(:,end);

if options.useOblique % correct for obliquity factor in this version
    nt_conv = length(wavelet);
    F = opFFTsym_conv(nt_conv);
    nf = size(F,1);
    Binv = opObliqInv(dt,nf,1);
    wavelet_invObliq = F' * Binv * F * wavelet(:);
end

    
    
wavelet(wavelet_window_end+1:end) = [];
wavelet(1:wavelet_window_start-1) = [];
if options.useOblique
    wavelet_invObliq(wavelet_window_end+1:end) = [];
    wavelet_invObliq(1:wavelet_window_start-1) = [];
end

% Determine start time and axis
lagtime = (wavelet_window_start-1-nt)*dt;
endtime = (wavelet_window_end-nt-1)*dt;
time_axis_wavelet = [lagtime:dt:endtime];

% Plot wavelet
figure('Name','Final estimated wavelet')
plot(time_axis_wavelet,wavelet)
xlabel('time (s)')
ylabel('amplitude')
title('Final estimated wavelet')

if options.useOblique
    figure('Name','Final estimated wavelet (corrected for obliquity factor/source-ghosting)')
    plot(time_axis_wavelet,wavelet_invObliq)
    xlabel('time (s)')
    ylabel('amplitude')
    title('Final estimated wavelet (corrected for obliquity factor/source-ghosting)')
end

%% Plot a shot record
figure('Name','Sample shot record at n_shots/2')

% Plot initial data
ah = subplot(1,4,1);
colormap('gray')
imagesc(preview_initialdata)
title('Initial data')
c_axis = get(ah,'clim');
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)

% Plot primaries
subplot(1,4,2)
colormap('gray')
imagesc(preview_primary, c_axis)
title('Est. primaries')
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)

% Plot IR
subplot(1,4,3)
colormap('gray')
imagesc(preview_primaryIR)
title('Est. Greens func (colorbar rescaled)')
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)

% Plot residual
subplot(1,4,4)
colormap('gray')
imagesc(preview_initialdata - preview_primary, c_axis)
title('Initial data minus est. primaries')
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)


%% Plot z-offset data
figure('Name','Zero-offset trace gather')

% Plot initial data
ah = subplot(1,3,1);
colormap('gray')
imagesc(zoffset_initialdata)
title('Initial data')
c_axis = get(ah,'clim');
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)

% Plot primaries
subplot(1,3,2)
colormap('gray')
imagesc(zoffset_primary, c_axis)
title('Est. primaries')
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)

% Plot residual
subplot(1,3,3)
colormap('gray')
imagesc(zoffset_initialdata - zoffset_primary, c_axis)
title('Initial data minus est. primaries')
xlabel('trace number')
ylabel('time (s)')
set(gca, 'YTick', time_axis)
set(gca, 'YTickLabel', time_axis_label)

