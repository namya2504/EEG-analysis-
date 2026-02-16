%% Preprocessing - bandpass filter and segment (Following Aljalal et al. 2022 Methodology)
% Following Aljalal et al. 2022 methodology
% - Bandpass filter: 0.5-32 Hz
% - 5th order Butterworth filter
% - Segment into 10-second epochs

clear all; close all; clc;

%% 1. Load
load('results/step1_loaded_data.mat');
fprintf('Loaded: %d channels, %d samples\n', size(data,1), size(data,2));

%% 2. Design bandpass filter
% 0.5-32 Hz bandpass, 5th order Butterworth (per paper)
fs = header.fs;
[b, a] = butter(5, [0.5 32]/(fs/2), 'bandpass');

%% 3. Apply filter to EEG channels
filtered_data = zeros(size(data));

for i = 1:length(eeg_idx)
    ch = eeg_idx(i);
    fprintf('Filtering %s\n', eeg_labels{i});
    
    % applying zero phase filter
    filtered_data(ch, :) = filtfilt(b, a, data(ch, :));
end

%% 4. Before/after comparison
figure('Position', [100, 100, 1400, 800]);

t_sec = 10;
n_samples = t_sec * fs;
t = (0:n_samples-1) / fs;
channel = eeg_idx(1);

% original
subplot(4,1,1)
plot(t, data(channel, 1:n_samples), 'b-', 'LineWidth', 0.8)
title('Original'); ylabel('μV'); grid on

% filtered
subplot(4,1,2)
plot(t, filtered_data(channel, 1:n_samples), 'r-', 'LineWidth', 0.8)
title('Filtered (0.5-32 Hz)'); ylabel('μV'); grid on
xlabel('Time (s)')

% PSD comparison
subplot(4,1,3)
[pxx_orig, f] = pwelch(data(channel, 1:n_samples), fs*4, fs*2, [], fs);
[pxx_filt, ~] = pwelch(filtered_data(channel, 1:n_samples), fs*4, fs*2, [], fs);
plot(f, 10*log10(pxx_orig), 'b-', 'LineWidth', 1.5); hold on
plot(f, 10*log10(pxx_filt), 'r-', 'LineWidth', 1.5);
xline([0.5 32], '--k', 'LineWidth', 1.5); hold off
xlim([0 50]); xlabel('Frequency (Hz)'); ylabel('Power (dB/Hz)');
title('Power Spectrum'); legend('Original', 'Filtered'); grid on

% filter response
subplot(4,1,4)
[h, w] = freqz(b, a, 1024, fs);
plot(w, 20*log10(abs(h)), 'k-', 'LineWidth', 2);
xline([0.5 32], '--r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Filter Response'); grid on; xlim([0 50]);

sgtitle(sprintf('Filtering: %s', eeg_labels{1}));

%% 5. Divide into 10 sec segments
seg_len = 10; % seconds
seg_samples = seg_len * fs;
n_segments = floor(size(filtered_data, 2) / seg_samples);

fprintf('\nSegmenting: %d samples -> %d epochs of %ds\n', ...
    size(filtered_data,2), n_segments, seg_len);

segments = zeros(length(eeg_idx), seg_samples, n_segments);

for seg = 1:n_segments
    start_idx = (seg-1) * seg_samples + 1;
    end_idx = seg * seg_samples;
    
    for i = 1:length(eeg_idx)
        ch = eeg_idx(i);
        segments(i, :, seg) = filtered_data(ch, start_idx:end_idx);
    end
end

fprintf('Final dims: %d channels × %d samples × %d segments\n', ...
    size(segments,1), size(segments,2), size(segments,3));

%% 6. Visualize sample segments
figure('Position', [100, 100, 1400, 600]);

seg_idx = randperm(n_segments, 3);
for i = 1:3
    subplot(3,1,i)
    t_seg = (0:seg_samples-1) / fs;
    plot(t_seg, squeeze(segments(1, :, seg_idx(i))), 'b-', 'LineWidth', 0.8);
    title(sprintf('Segment %d', seg_idx(i)));
    ylabel('μV'); grid on
    if i == 3, xlabel('Time (s)'); end
end

sgtitle('Sample Segments');

%% 7. Save
save('results/step2_preprocessed.mat', 'segments', 'eeg_idx', 'eeg_labels', ...
     'header', 'fs', 'seg_len', 'n_segments');

saveas(figure(1), 'figures/step2_filtering_comparison.png');
saveas(figure(2), 'figures/step2_sample_segments.png');

fprintf('\nPreprocessed data saved\n');
