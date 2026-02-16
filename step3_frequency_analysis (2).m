
%% STEP 3: Frequency Band Power Analysis

% Purpose: Extract power in different frequency bands to understand brain activity patterns during sleep
% pwelch fucntion used to estiamte PSD which is then integrated over 10 sec segments over each frequency band to find power int hat band
% Extract delta, theta, alpha, beta power from segments

clear all; close all; clc;

%% 1. Load preprocessed data from step 2
if ~isfile('results/step2_preprocessed.mat')
    error('Run step2_preprocess.m first');
end

load('results/step2_preprocessed.mat');
fprintf('Loaded %d segments, %d channels\n', n_segments, size(segments,1));

%% 2. Define bands
%These frequencies are based on existing literature on different stages of
%sleep
bands = struct();
bands.delta = [0.5, 4];
bands.theta = [4, 8];
bands.alpha = [8, 13];
bands.beta = [13, 30];
band_names = {'delta (Deep Sleep)', 'theta(Light sleep)', 'alpha(NREM, Drowsy)', 'beta(Active, REM)'};

%% 3. PSD parameters
% window size was chosen to be 4 to compromise between poor frequency resolution (due
% to shorter windows) or poor time resolution (due to longer windows)
win_size = fs * 4;  % 4 sec
overlap = win_size / 2;
nfft = win_size;

%% 4. Strorage of Power 
num_channels = size(segments, 1);
num_segments = size(segments, 3);
band_names = {'delta', 'theta', 'alpha', 'beta'};
num_bands = length(band_names);

band_power = struct();
for b = 1:num_bands
   
    band_power.(band_names{b}) = zeros(num_segments, num_channels);
end

%% 5. Extract power for all segments
fprintf('Extracting band power...\n');
tic;

for seg = 1:num_segments
    if mod(seg, 500) == 0
        fprintf('  %d/%d (%.0fs elapsed)\n', seg, num_segments, toc);
    end
    
    for ch = 1:num_channels
        signal = squeeze(segments(ch, :, seg));
        
        if any(isnan(signal)) || all(signal == 0)
            continue
        end
        
        [pxx, f] = pwelch(signal, hamming(win_size), overlap, nfft, fs);
        
        for b = 1:num_bands
            band_name = band_names{b};
            freq_range = bands.(band_name);
            freq_idx = (f >= freq_range(1)) & (f <= freq_range(2));
            band_power.(band_name)(seg, ch) = trapz(f(freq_idx), pxx(freq_idx));
        end
    end
end

fprintf('Done (%.1fs)\n\n', toc);

%% 6. Stats for primary channel
primary_channel = 1;
fprintf('Stats for %s:\n', eeg_labels{primary_channel});
fprintf('Band      Mean    Std     Min     Max     Median\n');
fprintf('-----------------------------------------------------\n');

stats = struct();
for b = 1:num_bands
    band_name = band_names{b};
    power_vals = band_power.(band_name)(:, primary_channel);
    
    stats.(band_name).mean = mean(power_vals);
    stats.(band_name).std = std(power_vals);
    stats.(band_name).min = min(power_vals);
    stats.(band_name).max = max(power_vals);
    stats.(band_name).median = median(power_vals);
    stats.(band_name).p25 = prctile(power_vals, 25);
    stats.(band_name).p75 = prctile(power_vals, 75);
    
    fprintf('%-8s  %.2f  %.2f  %.2f  %.2f  %.2f\n', ...
        band_name, stats.(band_name).mean, stats.(band_name).std, ...
        stats.(band_name).min, stats.(band_name).max, stats.(band_name).median);
end

%% 7. Time series plot
time_hours = (0:n_segments-1) * seg_len / 3600;

figure('Position', [100, 100, 1400, 900]);

for b = 1:num_bands
    subplot(4,1,b)
    band_name = band_names{b};
    power_vals = band_power.(band_name)(:, primary_channel);
    
    plot(time_hours, power_vals, 'LineWidth', 1.2);
    hold on
    yline(stats.(band_name).mean, '--r', 'LineWidth', 1.5);
    hold off
    
    ylabel('Power (μV²/Hz)');
    title(sprintf('%s (%.1f-%.0f Hz) - %s', upper(band_name), ...
        bands.(band_name)(1), bands.(band_name)(2), eeg_labels{primary_channel}));
    grid on
    ylim([0, max(power_vals) * 1.1]);
    
    if b == 4, xlabel('Time (hours)'); end
end

sgtitle('Band Power Over Time');


%% 8. Save
save('results/step3_band_powers.mat', 'band_power', 'stats', 'bands', ...
     'band_names', 'time_hours', 'primary_channel', 'eeg_labels', 'seg_len');

saveas(figure(1), 'figures/step3_power_timeseries.png');

fprintf('\nResults saved\n');
