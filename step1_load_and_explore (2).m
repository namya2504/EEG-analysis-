%% Load and check sleep EEG data
% This step will: 
% Load .edf file
% Look at what the data contains
% Make some basic plots 
% I obtained sleep PDG files from open access PhysioNet database.

clear all; close all; clc;

%% 1. Set up folders and paths
addpath('code');
data_folder = 'data';

%% 2. Get data File
files = dir(fullfile(data_folder, '*.edf'));
if isempty(files)
    error('No .edf files in data folder');
end

%% 3. Load data
filename = fullfile(data_folder, files(1).name);
fprintf('Loading %s...\n', files(1).name);
[data, header] = read_edf(filename);

fprintf('Duration: %.1f hrs, Fs=%d Hz, %d channels\n', ...
    size(data,2)/header.fs/3600, header.fs, header.channels);

%% 4. Find EEG channels
% look for Fpz-Cz or Pz-Oz type labels
eeg_idx = [];
eeg_labels = {};

for i = 1:header.channels
    lbl = upper(header.label{i});
    if contains(lbl, 'EEG') || contains(lbl, 'FPZ') || contains(lbl, 'PZ') || contains(lbl, 'CZ')
        eeg_idx = [eeg_idx; i];
        eeg_labels{end+1} = header.label{i};
    end
end

% fallback to first 2 channels if naming doesn't match
if isempty(eeg_idx)
    eeg_idx = [1; 2];
    eeg_labels = {header.label{1}, header.label{2}};
end

fprintf('Using channels: %s\n', strjoin(eeg_labels, ', '));

%% 5. Visualization - 30sec window
figure('Position', [100, 100, 1400, 900]);

t_show = 30;
n_show = min(t_show * header.fs, size(data,2));
time = (0:n_show-1) / header.fs;
signal = data(eeg_idx(1), 1:n_show);

% Plot 1: raw signal
subplot(3,1,1)
plot(time, signal, 'b-', 'LineWidth', 0.8)
xlabel('Time (s)'); ylabel('Amplitude (μV)');
title(sprintf('%s - first %ds', eeg_labels{1}, t_show));
grid on

% Plot 2: Power spectral density
subplot(3,1,2)
win = header.fs * 4;
overlap = header.fs * 2;
[pxx, f] = pwelch(signal, win, overlap, [], header.fs);
plot(f, 10*log10(pxx), 'r-', 'LineWidth', 1.5)
xlabel('Frequency (Hz)'); ylabel('Power (dB/Hz)');
title('Power Spectral Density');
xlim([0 50]); grid on

% Delta (0.5-4 Hz) = deep sleep
% Theta (4-8 Hz) = light sleep  
% Alpha (8-13 Hz) = drowsy
% Beta (13-30 Hz) = active/REM

hold on
xline(4, '--k', 'Delta|Theta');
xline(8, '--k', 'Theta|Alpha');
xline(13, '--k', 'Alpha|Beta');
xline(30, '--k', 'Beta');
hold off

% Plot 3: Time Frequency
subplot(3,1,3)
spectrogram(signal, win, overlap, [], header.fs, 'yaxis');
ylim([0 50]);
title('Time-Frequency');

sgtitle(files(1).name, 'Interpreter', 'none');

%% 6. Save
save('results/step1_loaded_data.mat', 'data', 'header', 'eeg_idx', 'eeg_labels');
saveas(gcf, 'figures/step1_initial_exploration.png');

fprintf('\nData saved to results/step1_loaded_data.mat\n');
