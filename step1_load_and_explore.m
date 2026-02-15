%% STEP 1: Load and Explore Sleep EEG Data
% Author: [Your Name]
% Date: [Date]
%
% This is my first time working with EEG data. I'm learning how to:
% - Load a .edf file (EEG recording format)
% - Look at what the data contains
% - Make some basic plots to see what brain waves look like
%
% I found the read_edf function online and modified it to work with
% sleep EEG files from the PhysioNet database.

clear all;  % Start fresh
close all;  % Close any open figures
clc;        % Clear the command window

fprintf('========================================\n');
fprintf('   LOADING SLEEP EEG DATA\n');
fprintf('========================================\n\n');

%% 1. Set up folders and paths
fprintf('Step 1: Setting up folders...\n');

% Make sure we're in the right directory
% (You might need to change this to where your project is)
% cd('PD_Sleep_Project')  % Uncomment and modify if needed

% Add the code folder so MATLAB can find our functions
addpath('code');

fprintf('  ✓ Folders set up\n\n');

%% 2. Find the EEG data file
fprintf('Step 2: Looking for data files...\n');

% My EEG files are in the 'data' folder
data_folder = 'data';

% List all .edf files (EEG file format)
files = dir(fullfile(data_folder, '*.edf'));

% Show what we found
fprintf('  Found %d EEG file(s):\n', length(files));
for i = 1:length(files)
    fprintf('    %d. %s\n', i, files(i).name);
end
fprintf('\n');

% Make sure we actually found something
if isempty(files)
    error('No EEG files found! Make sure .edf files are in the data/ folder');
end

%% 3. Load the first file
fprintf('Step 3: Loading the first file...\n');

% I'm going to analyze the first file
filename = fullfile(data_folder, files(1).name);
fprintf('  Loading: %s\n', files(1).name);
fprintf('  This might take 10-20 seconds...\n\n');

% Use the read_edf function to load the data
% This function reads the special .edf format that EEG uses
[data, header] = read_edf(filename);

fprintf('  ✓ File loaded!\n\n');

%% 4. Look at what we got
fprintf('Step 4: Examining the data...\n\n');

% The 'header' structure tells us about the recording
fprintf('  Recording Information:\n');
fprintf('  ---------------------\n');
fprintf('  Patient ID: %s\n', header.patient);
fprintf('  Date: %s\n', header.startdate);
fprintf('  Recording duration: %.1f hours\n', size(data,2)/header.fs/3600);
fprintf('  Sampling rate: %d Hz (samples per second)\n', header.fs);
fprintf('  Number of channels: %d\n\n', header.channels);

% Show what channels we have
fprintf('  Available channels:\n');
for i = 1:header.channels
    fprintf('    %d. %s\n', i, header.label{i});
end
fprintf('\n');

%% 5. Find the EEG channels
fprintf('Step 5: Finding EEG channels...\n');

% Sleep recordings have many channels (EEG, EOG, EMG, etc.)
% I only want the EEG channels for brain activity
% EEG channel names usually contain "EEG", "Fpz", "Pz", or "Cz"

eeg_idx = [];           % Will store channel numbers
eeg_labels = {};        % Will store channel names

for i = 1:header.channels
    % Convert to uppercase to make matching easier
    label_upper = upper(header.label{i});
    
    % Check if this looks like an EEG channel
    if contains(label_upper, 'EEG') || contains(label_upper, 'FPZ') || ...
       contains(label_upper, 'PZ') || contains(label_upper, 'CZ')
        
        % Found one! Add it to our list
        eeg_idx = [eeg_idx; i];
        eeg_labels{end+1} = header.label{i};
    end
end

fprintf('  Found %d EEG channel(s):\n', length(eeg_idx));
for i = 1:length(eeg_idx)
    fprintf('    - %s (channel #%d)\n', eeg_labels{i}, eeg_idx(i));
end
fprintf('\n');

% If we didn't find any EEG channels, just use the first 2 channels
if isempty(eeg_idx)
    fprintf('  Hmm, no channels named "EEG" found.\n');
    fprintf('  I''ll use the first 2 channels instead.\n\n');
    eeg_idx = [1; 2];
    eeg_labels = {header.label{1}, header.label{2}};
end

%% 6. Look at the actual brain wave data
fprintf('Step 6: Examining the signal...\n\n');

% Let's look at the first EEG channel
channel = eeg_idx(1);

% Get some basic statistics
signal_mean = mean(data(channel, :));
signal_std = std(data(channel, :));
signal_min = min(data(channel, :));
signal_max = max(data(channel, :));

fprintf('  Channel: %s\n', eeg_labels{1});
fprintf('  Average value: %.2f μV (microvolts)\n', signal_mean);
fprintf('  Variability (std): %.2f μV\n', signal_std);
fprintf('  Range: %.2f to %.2f μV\n\n', signal_min, signal_max);

% These numbers tell me the signal looks reasonable
% EEG is typically ±100 μV or so

%% 7. Create visualizations
fprintf('Step 7: Creating plots...\n');
fprintf('  Making 3 different views of the data...\n\n');

% Create a new figure window
figure('Position', [100, 100, 1400, 900], 'Name', 'Sleep EEG Exploration');

% I'll show 30 seconds of data (easier to see patterns than hours)
time_to_show = 30;  % seconds
samples_to_show = time_to_show * header.fs;

% Make sure we don't try to show more data than we have
if samples_to_show > size(data, 2)
    samples_to_show = size(data, 2);
    time_to_show = samples_to_show / header.fs;
end

% Create a time vector (for the x-axis)
time = (0:samples_to_show-1) / header.fs;

% Get the signal for these 30 seconds
signal = data(channel, 1:samples_to_show);

%--- PLOT 1: Raw brain waves ---
subplot(3, 1, 1)
plot(time, signal, 'b-', 'LineWidth', 0.8)
xlabel('Time (seconds)')
ylabel('Amplitude (μV)')
title(sprintf('Raw EEG Signal: %s (first %d seconds)', eeg_labels{1}, time_to_show), ...
      'FontWeight', 'bold')
grid on

% The wiggly line shows brain electrical activity!

%--- PLOT 2: Frequency content ---
subplot(3, 1, 2)

% Use pwelch to see which frequencies are in the signal
% This is a built-in MATLAB function for frequency analysis
window_length = header.fs * 4;  % 4 second windows
overlap = header.fs * 2;         % 50% overlap (standard)

% Calculate power spectral density
[power, frequency] = pwelch(signal, window_length, overlap, [], header.fs);

% Plot it
plot(frequency, 10*log10(power), 'r-', 'LineWidth', 1.5)
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
title('Power Spectrum - Which frequencies are present?', 'FontWeight', 'bold')
grid on
xlim([0 50])  % Only show 0-50 Hz (where brain waves are)

% Add lines to show different brain wave bands
% I learned these from sleep research papers:
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

%--- PLOT 3: Time-frequency view ---
subplot(3, 1, 3)

% Spectrogram shows how frequencies change over time
% (combines time view + frequency view)
spectrogram(signal, window_length, overlap, [], header.fs, 'yaxis');
title('Spectrogram - Frequencies over time', 'FontWeight', 'bold')
ylim([0 50])  % Only show 0-50 Hz

% Add overall title
sgtitle(sprintf('Sleep EEG Data: %s', files(1).name), ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('  ✓ Plots created!\n\n');

%% 8. Save the data for next steps
fprintf('Step 8: Saving data...\n');

% Save the important stuff so I don't have to reload the big file
save('results/step1_loaded_data.mat', 'data', 'header', 'eeg_idx', 'eeg_labels');

fprintf('  ✓ Saved to: results/step1_loaded_data.mat\n\n');

% Also save the figure
saveas(gcf, 'figures/step1_initial_exploration.png');
fprintf('  ✓ Figure saved to: figures/step1_initial_exploration.png\n\n');

%% 9. Summary
fprintf('========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n\n');

fprintf('What I found:\n');
fprintf('  - Loaded %.1f hours of sleep data\n', size(data,2)/header.fs/3600);
fprintf('  - Found %d EEG channel(s)\n', length(eeg_idx));
fprintf('  - Signal looks normal (range: %.1f to %.1f μV)\n', signal_min, signal_max);
fprintf('  - Can see brain wave patterns in the plots\n\n');

fprintf('What the plots show:\n');
fprintf('  1. Top: Wiggly line = brain electrical activity\n');
fprintf('  2. Middle: Which frequencies are strongest\n');
fprintf('  3. Bottom: How frequencies change over time\n\n');

fprintf('Next step:\n');
fprintf('  Run step2_preprocess.m to clean up the signal\n\n');

fprintf('========================================\n');
fprintf('✓ STEP 1 COMPLETE!\n');
fprintf('========================================\n\n');