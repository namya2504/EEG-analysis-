%% Step 4 - Looking for patterns in the sleep data
% Purpose: Find unusual segments and assess sleep quality using differences
% between early and late sleep

clear all; close all; clc;

fprintf('Starting analysis...\n\n');

%% 1. Load data calculated in step 3
if ~isfile('results/step3_band_powers.mat')
    error('Run Step 3');
end

load('results/step3_band_powers.mat');

fprintf('Loaded data for %s\n', eeg_labels{primary_channel});
fprintf('Total recording time: %.1f hours\n\n', max(time_hours));

%% 2. Find unusual segments

fprintf('Looking for unusual segments...\n');

% Get the power values for each band
delta_vals = band_power.delta(:, primary_channel);
beta_vals = band_power.beta(:, primary_channel);
theta_vals = band_power.theta(:, primary_channel);
alpha_vals = band_power.alpha(:, primary_channel);

% Find unsual values (top 5%)
high_delta = find(delta_vals > prctile(delta_vals, 95));
high_beta = find(beta_vals > prctile(beta_vals, 95));

% Unusual segments have BOTH high at the same time
unusual_segments = intersect(high_delta, high_beta);

fprintf('Found %d segments with high delta\n', length(high_delta));
fprintf('Found %d segments with high beta\n', length(high_beta));
fprintf('Found %d segments with BOTH high \n\n', length(unusual_segments));

%% 3. Try to plot these unusual segments

if ~isempty(unusual_segments)
    fprintf('Plotting unusual segments...\n');
    
    step2_data = load('results/step2_preprocessed.mat');
    segs = step2_data.segments;
    fs = step2_data.fs;
    

    if isfield(step2_data, 'segment_length')
        seg_len = step2_data.segment_length;
    else
        seg_len = 10;  
    end
    
    fprintf('Loaded %d segments\n', size(segs, 3));
    
 
    unusual_segments = unusual_segments(unusual_segments <= size(segs, 3));
    
    if ~isempty(unusual_segments)
        % Create a figure to show them
        figure('Position', [100, 100, 1400, 900]);
        
        % Show up to 3 examples
        num_show = min(3, length(unusual_segments));
        
        for i = 1:num_show
            idx = unusual_segments(i);
            
            % Get the signal for this segment
            signal = squeeze(segs(primary_channel, :, idx));
            time = (0:length(signal)-1) / fs;
            
            % Left plot: actual brain waves
            subplot(num_show, 2, (i-1)*2 + 1)
            plot(time, signal, 'b-');
            xlabel('Time (sec)');
            ylabel('Amplitude (uV)');
            title(sprintf('Segment %d (%.1f hours)', idx, time_hours(idx)));
            grid on;
            
            % Right plot: compare to average
            subplot(num_show, 2, (i-1)*2 + 2)
            
            % This segment's powers
            seg_powers = [delta_vals(idx), theta_vals(idx), ...
                         alpha_vals(idx), beta_vals(idx)];
            
            % Average powers across all segments
            avg_powers = [mean(delta_vals), mean(theta_vals), ...
                         mean(alpha_vals), mean(beta_vals)];
            
            % Plotting
            data_to_plot = [avg_powers; seg_powers];
            bar(data_to_plot');
            
            set(gca, 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta'});
            ylabel('Power');
            title('This segment vs average');
            legend({'Average', 'This one'});
            grid on;
        end
        
        sgtitle('Unusual Segments (high delta AND high beta)');
        fprintf('Created plot\n\n');
    end
else
    fprintf('No unusual segments found\n\n');
end

%% 4. Compare first half of night vs second half
% Disruptions in this can indicate poor quality of sleep

fprintf('Comparing early vs late sleep...\n');

% Split the data in half
mid = floor(length(time_hours) / 2);

first_half = band_power.delta(1:mid, primary_channel);
second_half = band_power.delta(mid+1:end, primary_channel);

mean1 = mean(first_half);
mean2 = mean(second_half);

fprintf('First half mean: %.2f\n', mean1);
fprintf('Second half mean: %.2f\n', mean2);

diff = mean1 - mean2;
fprintf('Difference: %.2f\n', diff);

if diff > 0
    fprintf('-> More delta in first half (normal pattern)\n\n');
else
    fprintf('-> More delta in second half (unusual?)\n\n');
end

%% 5. Create plots

fprintf('Making comparison plots...\n');

% Plot 1: bar chart comparing first vs second half
figure;

subplot(1,2,1)
bar([mean1, mean2]);
set(gca, 'XTickLabel', {'First half', 'Second half'});
ylabel('Mean delta power');
title('First vs second half of night');
grid on;

% boxplot 
subplot(1,2,2)
boxplot([first_half; second_half], [ones(size(first_half)); 2*ones(size(second_half))]);
set(gca, 'XTickLabel', {'First', 'Second'});
ylabel('Delta power');
title('Distribution');
grid on;

% Plot 2: all the frequency bands 
figure;
all_data = [delta_vals, theta_vals, alpha_vals, beta_vals];
boxplot(all_data);
set(gca, 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta'});
ylabel('Power');
title('All frequency bands');
grid on;

fprintf('Done with plots\n\n');

%% 6. Saving Results

fprintf('Saving results...\n');


results = struct();
results.period_comparison.mean_first = mean1;
results.period_comparison.mean_second = mean2;
results.extreme_segments.high_delta = high_delta;
results.extreme_segments.high_beta = high_beta;
results.extreme_segments.unusual = unusual_segments;

save('results/step4_exploration.mat', 'results');

% Save the figures
saveas(figure(1), 'figures/step4_period_comparison.png');
saveas(figure(2), 'figures/step4_band_distributions.png');

fprintf('Saved Results\n');
fprintf('Step 4 complete\n');