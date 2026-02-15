%% STEP 3: Frequency Band Power Analysis

% Purpose: Extract power in different frequency bands to understand
%          brain activity patterns during sleep
%
% BACKGROUND:
% Different brain wave frequencies correspond to different brain states:
%   - Delta (0.5-4 Hz): Deep sleep, slow-wave activity
%   - Theta (4-8 Hz): Light sleep, drowsiness, REM transitions
%   - Alpha (8-13 Hz): Relaxed wakefulness, drowsiness
%   - Beta (13-30 Hz): Active thinking, REM sleep
%
% METHODS:
% I'm using MATLAB's pwelch function to estimate the power spectral density,
% then integrating over each frequency band to get total power in that band
% using 10 second segments.

clear all; close all; clc;

fprintf('========================================\n');
fprintf('  FREQUENCY BAND POWER ANALYSIS\n');
fprintf('========================================\n\n');

%% 1. Load preprocessed data from Step 2
fprintf('STEP 1: Loading preprocessed data...\n');

% Check if file exists
if ~isfile('../results/step2_preprocessed.mat')
    error('Step 2 data not found. Please run step2_preprocess.m first.');
end

load('../results/step2_preprocessed.mat');

fprintf('✓ Loaded preprocessed data:\n');
fprintf('  - Segments: %d\n', num_segments);
fprintf('  - Channels: %d\n', size(segments, 1));
fprintf('  - Samples per segment: %d\n', size(segments, 2));
fprintf('  - Sampling rate: %d Hz\n', fs);
fprintf('  - Segment duration: %d seconds\n\n', segment_length);

%% 2. Define frequency bands based on neuroscience literature
fprintf('STEP 2: Defining frequency bands...\n');

% These frequency ranges are standard in sleep EEG research
bands = struct();
bands.delta = [0.5, 4];    % Deep sleep
bands.theta = [4, 8];      % Light sleep, REM
bands.alpha = [8, 13];     % Drowsiness, relaxed
bands.beta = [13, 30];     % Active, REM

fprintf('  Frequency bands defined:\n');
fprintf('    Delta: %.1f-%.0f Hz (deep sleep)\n', bands.delta(1), bands.delta(2));
fprintf('    Theta: %.0f-%.0f Hz (light sleep, REM)\n', bands.theta(1), bands.theta(2));
fprintf('    Alpha: %.0f-%.0f Hz (drowsiness)\n', bands.alpha(1), bands.alpha(2));
fprintf('    Beta: %.0f-%.0f Hz (active, REM)\n\n', bands.beta(1), bands.beta(2));

%% 3. Set up parameters for power spectral density estimation
fprintf('STEP 3: Setting up PSD parameters...\n');

% Window size for Welch's method
% I chose 4 seconds because:
% - Shorter windows = poor frequency resolution (can't separate bands)
% - Longer windows = poor time resolution (can't track changes)
% - 4 seconds is a good compromise for sleep EEG
window_size = fs * 4;  % 4 seconds = 400 samples at 100 Hz

% Overlap: 50% is standard (balances smoothing and independence)
overlap_size = window_size / 2;

% FFT points: same as window for computational efficiency
nfft = window_size;

fprintf('  PSD parameters:\n');
fprintf('    Window size: %d samples (%.1f seconds)\n', window_size, window_size/fs);
fprintf('    Overlap: %d samples (%.1f seconds)\n', overlap_size, overlap_size/fs);
fprintf('    Frequency resolution: %.2f Hz\n\n', fs/nfft);

%% 4. Initialize storage arrays
fprintf('STEP 4: Initializing storage arrays...\n');

% We'll store power for each band, each channel, each segment
% Dimensions: [num_segments × num_channels × num_bands]
num_channels = size(segments, 1);
band_names = {'delta', 'theta', 'alpha', 'beta'};
num_bands = length(band_names);

% Create structure to store band powers
band_power = struct();
for b = 1:num_bands
    % Initialize matrix: segments × channels
    band_power.(band_names{b}) = zeros(num_segments, num_channels);
end

fprintf('✓ Storage arrays initialized\n');
fprintf('  Will extract %d bands × %d channels × %d segments = %d values\n\n', ...
        num_bands, num_channels, num_segments, num_bands*num_channels*num_segments);

%% 5. Extract band power for all segments
fprintf('STEP 5: Extracting band power for all segments...\n');
fprintf('  This will take 1-2 minutes...\n\n');

% Start timer to track progress
tic;

% Loop through all segments
for seg = 1:num_segments
    
    % Progress update every 500 segments
    if mod(seg, 500) == 0
        elapsed = toc;
        percent_done = 100 * seg / num_segments;
        est_total = elapsed / (seg / num_segments);
        est_remaining = est_total - elapsed;
        
        fprintf('  Progress: %d/%d segments (%.1f%%) - %.0f sec elapsed, ~%.0f sec remaining\n', ...
                seg, num_segments, percent_done, elapsed, est_remaining);
    end
    
    % Loop through each channel
    for ch = 1:num_channels
        
        % Extract this segment
        signal = squeeze(segments(ch, :, seg));
        
        % Check for any issues with the signal
        if any(isnan(signal)) || all(signal == 0)
            warning('Segment %d, Channel %d has NaN or all zeros - skipping', seg, ch);
            continue;
        end
        
        % Compute power spectral density using Welch's method
        % This is more robust than raw FFT because it averages multiple windows
        [pxx, f] = pwelch(signal, hamming(window_size), overlap_size, nfft, fs);
        
        % Extract power in each frequency band
        % Method: Integrate PSD over the frequency range
        for b = 1:num_bands
            band_name = band_names{b};
            freq_range = bands.(band_name);
            
            % Find frequency indices that fall within this band
            freq_idx = (f >= freq_range(1)) & (f <= freq_range(2));
            
            % Integrate PSD over frequency range (trapezoid rule)
            % This gives total power in μV²/Hz
            band_power.(band_name)(seg, ch) = trapz(f(freq_idx), pxx(freq_idx));
        end
    end
end

total_time = toc;

fprintf('\n✓ Band power extraction complete!\n');
fprintf('  Total time: %.1f seconds (%.2f segments/sec)\n\n', ...
        total_time, num_segments/total_time);

%% 6. Compute basic statistics for each band
fprintf('STEP 6: Computing statistics...\n\n');

% We'll focus on the first channel (Fpz-Cz) for primary analysis
% Channel 1 is frontal-central, which is implicated in PD research
primary_channel = 1;

fprintf('  Primary analysis channel: %s\n\n', eeg_labels{primary_channel});

% Create statistics structure
stats = struct();

fprintf('  Band      Mean      Std       Min       Max      25th%%     50th%%     75th%%\n');
fprintf('  ------------------------------------------------------------------------\n');

for b = 1:num_bands
    band_name = band_names{b};
    
    % Get power values for primary channel
    power_values = band_power.(band_name)(:, primary_channel);
    
    % Compute statistics
    stats.(band_name).mean = mean(power_values);
    stats.(band_name).std = std(power_values);
    stats.(band_name).min = min(power_values);
    stats.(band_name).max = max(power_values);
    stats.(band_name).median = median(power_values);
    stats.(band_name).p25 = prctile(power_values, 25);
    stats.(band_name).p75 = prctile(power_values, 75);
    
    % Display
    fprintf('  %-8s  %7.2f   %7.2f   %7.2f   %7.2f   %7.2f   %7.2f   %7.2f\n', ...
            band_name, ...
            stats.(band_name).mean, ...
            stats.(band_name).std, ...
            stats.(band_name).min, ...
            stats.(band_name).max, ...
            stats.(band_name).p25, ...
            stats.(band_name).median, ...
            stats.(band_name).p75);
end

fprintf('\n');

%% 7. Create time series visualization
fprintf('STEP 7: Creating visualizations...\n');

% Create time vector (in hours)
time_hours = (0:num_segments-1) * segment_length / 3600;

% FIGURE 1: Band power over time
figure('Position', [100, 100, 1400, 900], 'Name', 'Band Power Time Series');

for b = 1:num_bands
    subplot(4, 1, b)
    
    band_name = band_names{b};
    power_values = band_power.(band_name)(:, primary_channel);
    
    % Plot time series
    plot(time_hours, power_values, 'LineWidth', 1.2)
    
    % Add mean line
    hold on
    yline(stats.(band_name).mean, '--r', 'LineWidth', 1.5, ...
          'Label', sprintf('Mean=%.1f', stats.(band_name).mean));
    hold off
    
    % Labels and formatting
    ylabel('Power (μV²/Hz)', 'FontSize', 11)
    title(sprintf('%s Band (%.1f-%.0f Hz) - Channel: %s', ...
                  upper(band_name), bands.(band_name)(1), bands.(band_name)(2), ...
                  eeg_labels{primary_channel}), ...
          'FontSize', 12, 'FontWeight', 'bold')
    grid on
    
    if b == 4
        xlabel('Time (hours)', 'FontSize', 11)
    end
    
    % Set consistent y-axis limits for better comparison
    ylim([0, max(power_values) * 1.1])
end

sgtitle('Frequency Band Power Evolution Over Sleep Period', ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('  ✓ Time series plot created\n');

%% 8. Create comparison bar chart
fprintf('  Creating comparison plots...\n');

% FIGURE 2: Comparison across bands
figure('Position', [100, 100, 1200, 600], 'Name', 'Band Power Comparison');

subplot(1, 2, 1)
% Bar chart with error bars
means = [stats.delta.mean, stats.theta.mean, stats.alpha.mean, stats.beta.mean];
stds = [stats.delta.std, stats.theta.std, stats.alpha.std, stats.beta.std];

bar_handle = bar(means, 'FaceColor', [0.3 0.5 0.8]);
hold on
errorbar(1:4, means, stds, 'k.', 'LineWidth', 1.5)
hold off

set(gca, 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta'})
ylabel('Mean Power (μV²/Hz)', 'FontSize', 11)
title('Average Power by Frequency Band', 'FontSize', 12, 'FontWeight', 'bold')
grid on

subplot(1, 2, 2)
% Box plots showing distribution
power_matrix = [band_power.delta(:, primary_channel), ...
                band_power.theta(:, primary_channel), ...
                band_power.alpha(:, primary_channel), ...
                band_power.beta(:, primary_channel)];

boxplot(power_matrix, 'Labels', {'Delta', 'Theta', 'Alpha', 'Beta'})
ylabel('Power (μV²/Hz)', 'FontSize', 11)
title('Distribution of Power Values', 'FontSize', 12, 'FontWeight', 'bold')
grid on

sgtitle('Frequency Band Power Comparison', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('  ✓ Comparison plots created\n');

%% 9. Save results
fprintf('\nSTEP 8: Saving results...\n');

% Save band power data and statistics
save('../results/step3_band_powers.mat', 'band_power', 'stats', 'bands', ...
     'band_names', 'time_hours', 'primary_channel', 'eeg_labels');

fprintf('  ✓ Data saved to: ../results/step3_band_powers.mat\n');

% Save figures
saveas(figure(1), '../figures/step3_power_timeseries.png');
saveas(figure(2), '../figures/step3_band_comparison.png');

fprintf('  ✓ Figures saved to figures/ folder\n');

%% 10. Summary of key observations
fprintf('\n========================================\n');
fprintf('KEY OBSERVATIONS:\n');
fprintf('========================================\n\n');

% Find when delta power is highest
[max_delta, max_delta_idx] = max(band_power.delta(:, primary_channel));
max_delta_time = time_hours(max_delta_idx);

fprintf('1. Delta power (deep sleep indicator):\n');
fprintf('   - Mean: %.2f μV²/Hz\n', stats.delta.mean);
fprintf('   - Peak at: %.1f hours (power = %.2f)\n', max_delta_time, max_delta);
fprintf('   - Interpretation: Highest delta typically occurs in first sleep cycles\n\n');

% Compare first vs second half (preliminary observation)
midpoint = floor(num_segments / 2);
delta_first_half = mean(band_power.delta(1:midpoint, primary_channel));
delta_second_half = mean(band_power.delta(midpoint+1:end, primary_channel));

fprintf('2. Delta power comparison:\n');
fprintf('   - First half of night: %.2f μV²/Hz\n', delta_first_half);
fprintf('   - Second half of night: %.2f μV²/Hz\n', delta_second_half);
fprintf('   - Difference: %.2f μV²/Hz (%.1f%% change)\n', ...
        delta_first_half - delta_second_half, ...
        100 * (delta_first_half - delta_second_half) / delta_first_half);
if delta_first_half > delta_second_half
    fprintf('   - Observation: More deep sleep in first half (expected pattern)\n\n');
else
    fprintf('   - Observation: Unusual pattern - would need to investigate\n\n');
end

% Beta observations
fprintf('3. Beta power (REM/arousal indicator):\n');
fprintf('   - Mean: %.2f μV²/Hz\n', stats.beta.mean);
fprintf('   - Variability (CV): %.2f\n', stats.beta.std / stats.beta.mean);
fprintf('   - Note: Visual inspection of time series may reveal periodic patterns\n');
fprintf('         (look for ~90-minute cycles suggesting REM periods)\n\n');

% Overall variability
fprintf('4. Relative variability across bands:\n');
for b = 1:num_bands
    band_name = band_names{b};
    cv = stats.(band_name).std / stats.(band_name).mean;
    fprintf('   - %s: CV = %.2f\n', band_name, cv);
end
fprintf('   - Higher CV suggests more dynamic changes in that frequency range\n\n');

fprintf('========================================\n');
fprintf('✓✓✓ STEP 3 COMPLETE! ✓✓✓\n');
fprintf('========================================\n\n');

fprintf('NEXT STEPS:\n');
fprintf('  - Run step4_exploratory_patterns.m for statistical analysis\n');
fprintf('  - Examine time series plots for periodic patterns\n');
fprintf('  - Compare patterns to known sleep physiology\n\n');





            