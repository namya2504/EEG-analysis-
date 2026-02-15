%% STEP 2: Preprocessing (Following Aljalal et al. 2022 Methodology)
% This applies the EXACT preprocessing pipeline from the paper:
% - Bandpass filter: 0.5-32 Hz
% - 5th order Butterworth filter
% - Segment into 10-second epochs

clear all; close all; clc;

fprintf('========================================\n');
fprintf('   SLEEP EEG PREPROCESSING - STEP 2\n');
fprintf('   (Following Aljalal et al. 2022)\n');
fprintf('========================================\n\n');

%% 1. Load data from Step 1
fprintf('STEP 1: Loading data from previous step...\n');
load('../results/step1_loaded_data.mat');
fprintf('✓ Loaded: %d channels, %d samples\n', size(data,1), size(data,2));

%% 2. Design Bandpass Filter
fprintf('\nSTEP 2: Designing bandpass filter...\n');
fprintf('  Frequency range: 0.5-32 Hz (as per Aljalal et al.)\n');
fprintf('  Filter order: 5 (Butterworth)\n');

fs = header.fs;  % Sampling rate (100 Hz)
low_freq = 0.5;  % Low cutoff (removes DC drift)
high_freq = 32;  % High cutoff (removes high-freq noise)
filter_order = 5;

% Design the filter
[b, a] = butter(filter_order, [low_freq high_freq]/(fs/2), 'bandpass');

fprintf('✓ Filter designed!\n');

%% 3. Apply Filter to EEG Channels
fprintf('\nSTEP 3: Applying filter to EEG channels...\n');

filtered_data = zeros(size(data));

for i = 1:length(eeg_idx)
    ch = eeg_idx(i);
    fprintf('  Filtering channel %d/%d: %s\n', i, length(eeg_idx), eeg_labels{i});
    
    % Apply zero-phase filter (filtfilt = forward + backward, no delay)
    filtered_data(ch, :) = filtfilt(b, a, data(ch, :));
end

fprintf('✓ Filtering complete!\n');

%% 4. Visualize: Before vs After Filtering
fprintf('\nSTEP 4: Creating before/after comparison...\n');

figure('Position', [100, 100, 1400, 800], 'Name', 'Filtering Results');

% Take 10 seconds for comparison
time_window = 10;
samples = time_window * fs;
t = (0:samples-1) / fs;

channel = eeg_idx(1);  % First EEG channel

% Plot original
subplot(4, 1, 1)
plot(t, data(channel, 1:samples), 'b-', 'LineWidth', 0.8)
title('ORIGINAL (Raw Signal)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Amplitude (µV)')
grid on

% Plot filtered
subplot(4, 1, 2)
plot(t, filtered_data(channel, 1:samples), 'r-', 'LineWidth', 0.8)
title('FILTERED (0.5-32 Hz Bandpass)', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Time (seconds)')
ylabel('Amplitude (µV)')
grid on

% Plot power spectra comparison
subplot(4, 1, 3)
[pxx_orig, f] = pwelch(data(channel, 1:samples), fs*4, fs*2, [], fs);
[pxx_filt, ~] = pwelch(filtered_data(channel, 1:samples), fs*4, fs*2, [], fs);

plot(f, 10*log10(pxx_orig), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(f, 10*log10(pxx_filt), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered')
xline(0.5, '--k', 'Low cutoff', 'LineWidth', 1.5)
xline(32, '--k', 'High cutoff', 'LineWidth', 1.5)
hold off
xlim([0 50])
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
title('POWER SPECTRUM: Before vs After', 'FontSize', 12, 'FontWeight', 'bold')
legend('Location', 'northeast')
grid on

% Plot frequency response of filter
subplot(4, 1, 4)
[h, w] = freqz(b, a, 1024, fs);
plot(w, 20*log10(abs(h)), 'k-', 'LineWidth', 2)
xline(0.5, '--r', 'LineWidth', 1.5)
xline(32, '--r', 'LineWidth', 1.5)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('FILTER FREQUENCY RESPONSE', 'FontSize', 12, 'FontWeight', 'bold')
grid on
xlim([0 50])

sgtitle(sprintf('Preprocessing: Bandpass Filtering (Channel: %s)', eeg_labels{1}), ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('✓ Comparison plots created!\n');

%% 5. Segment into 10-second Epochs (Aljalal method)
fprintf('\nSTEP 5: Segmenting into epochs...\n');
fprintf('  Epoch length: 10 seconds (as per Aljalal et al.)\n');

segment_length = 10;  % seconds
segment_samples = segment_length * fs;

% Calculate how many complete segments we can make
total_samples = size(filtered_data, 2);
num_segments = floor(total_samples / segment_samples);

fprintf('  Total samples: %d\n', total_samples);
fprintf('  Samples per segment: %d\n', segment_samples);
fprintf('  Number of segments: %d\n', num_segments);

% Create 3D array: [channels × samples × segments]
segments = zeros(length(eeg_idx), segment_samples, num_segments);

for seg = 1:num_segments
    if mod(seg, 100) == 0
        fprintf('  Processing segment %d/%d\n', seg, num_segments);
    end
    
    start_idx = (seg-1) * segment_samples + 1;
    end_idx = seg * segment_samples;
    
    for i = 1:length(eeg_idx)
        ch = eeg_idx(i);
        segments(i, :, seg) = filtered_data(ch, start_idx:end_idx);
    end
end

fprintf('✓ Segmentation complete!\n');
fprintf('  Final dimensions: %d channels × %d samples × %d segments\n', ...
        size(segments, 1), size(segments, 2), size(segments, 3));

%% 6. Visualize Some Segments
fprintf('\nSTEP 6: Visualizing sample segments...\n');

figure('Position', [100, 100, 1400, 600], 'Name', 'Sample Segments');

% Show 3 random segments
seg_indices = randperm(num_segments, 3);

for i = 1:3
    subplot(3, 1, i)
    t_seg = (0:segment_samples-1) / fs;
    plot(t_seg, squeeze(segments(1, :, seg_indices(i))), 'b-', 'LineWidth', 0.8)
    title(sprintf('Segment #%d (10 seconds)', seg_indices(i)), ...
          'FontSize', 11, 'FontWeight', 'bold')
    ylabel('Amplitude (µV)')
    grid on
    
    if i == 3
        xlabel('Time (seconds)')
    end
end

sgtitle('Random Preprocessed Segments', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('✓ Segment visualization created!\n');

%% 7. Save Preprocessed Data
fprintf('\nSTEP 7: Saving preprocessed data...\n');

save('../results/step2_preprocessed.mat', 'segments', 'eeg_idx', 'eeg_labels', ...
     'header', 'fs', 'segment_length', 'num_segments');

fprintf('✓ Saved to: results/step2_preprocessed.mat\n');

% Save figures
saveas(figure(1), '../figures/step2_filtering_comparison.png');
saveas(figure(2), '../figures/step2_sample_segments.png');

fprintf('✓ Figures saved to figures/ folder\n');

%% 8. Summary Statistics
fprintf('\n========================================\n');
fprintf('PREPROCESSING SUMMARY:\n');
fprintf('========================================\n');
fprintf('Input:\n');
fprintf('  - Raw data: %d channels × %d samples\n', size(data,1), size(data,2));
fprintf('  - Duration: %.1f hours\n', size(data,2)/fs/3600);
fprintf('\nProcessing:\n');
fprintf('  - Filter: 0.5-32 Hz bandpass (5th order Butterworth)\n');
fprintf('  - Segmentation: %d-second epochs\n', segment_length);
fprintf('\nOutput:\n');
fprintf('  - Preprocessed segments: %d channels × %d samples × %d segments\n', ...
          size(segments,1), size(segments,2), size(segments,3));
fprintf('  - Total epochs: %d (%.1f minutes of data)\n', ...
          num_segments, num_segments*segment_length/60);

fprintf('\n========================================\n');
fprintf('STEP 2 COMPLETE\n');
fprintf('========================================\n');
fprintf('\nData is  for Step 3\n');

fprintf('\nNext step: Band Power Analysis \n');