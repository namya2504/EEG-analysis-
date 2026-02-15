%% STEP 4: Unusual Segment Detection
% Purpose: Identify and visualize segments with unusual/disturbed patterns
clear all; close all; clc;

fprintf('========================================\n');
fprintf('  UNUSUAL SEGMENT DETECTION\n');
fprintf('========================================\n\n');

%% 1. Load band power data from Step 3
fprintf('STEP 1: Loading band power data...\n');

if ~isfile('../results/step3_band_powers.mat')
    error('Step 3 data not found. Please run step3_frequency_analysis.m first.');
end

load('../results/step3_band_powers.mat');

fprintf('✓ Data loaded successfully\n');
fprintf('  - Analyzed channel: %s\n', eeg_labels{primary_channel});
fprintf('  - Number of segments: %d\n', length(time_hours));
fprintf('  - Duration: %.1f hours\n\n', max(time_hours));

%% 2. Identify Unusual Segments
fprintf('========================================\n');
fprintf('IDENTIFYING UNUSUAL SEGMENTS\n');
fprintf('========================================\n\n');

fprintf('Looking for segments with unusual patterns...\n\n');

% Get power values for primary channel
delta_values = band_power.delta(:, primary_channel);
theta_values = band_power.theta(:, primary_channel);
alpha_values = band_power.alpha(:, primary_channel);
beta_values = band_power.beta(:, primary_channel);

%--- Type 1: Extremely Low Delta (Insufficient Deep Sleep) ---
fprintf('Type 1: Low Delta Segments\n');
delta_p05 = prctile(delta_values, 5);
low_delta_idx = find(delta_values < delta_p05);

fprintf('  Threshold: %.2f μV²/Hz (5th percentile)\n', delta_p05);
fprintf('  Found: %d segments (%.1f%% of total)\n', ...
        length(low_delta_idx), 100*length(low_delta_idx)/length(delta_values));

if ~isempty(low_delta_idx)
    fprintf('  Times: ');
    if length(low_delta_idx) <= 5
        fprintf('%.1fh ', time_hours(low_delta_idx));
    else
        fprintf('%.1fh, %.1fh, %.1fh ... (showing first 3)', ...
                time_hours(low_delta_idx(1)), time_hours(low_delta_idx(2)), time_hours(low_delta_idx(3)));
    end
    fprintf('\n');
    fprintf('  Example (segment #%d at %.1fh):\n', low_delta_idx(1), time_hours(low_delta_idx(1)));
    fprintf('    Delta=%.2f, Theta=%.2f, Alpha=%.2f, Beta=%.2f\n', ...
            delta_values(low_delta_idx(1)), theta_values(low_delta_idx(1)), ...
            alpha_values(low_delta_idx(1)), beta_values(low_delta_idx(1)));
end
fprintf('\n');

%--- Type 2: Extremely High Beta (Excessive Arousal) ---
fprintf('Type 2: High Beta Segments\n');
beta_p95 = prctile(beta_values, 95);
high_beta_idx = find(beta_values > beta_p95);

fprintf('  Threshold: %.2f μV²/Hz (95th percentile)\n', beta_p95);
fprintf('  Found: %d segments (%.1f%% of total)\n', ...
        length(high_beta_idx), 100*length(high_beta_idx)/length(beta_values));

if ~isempty(high_beta_idx)
    fprintf('  Times: ');
    if length(high_beta_idx) <= 5
        fprintf('%.1fh ', time_hours(high_beta_idx));
    else
        fprintf('%.1fh, %.1fh, %.1fh ... (showing first 3)', ...
                time_hours(high_beta_idx(1)), time_hours(high_beta_idx(2)), time_hours(high_beta_idx(3)));
    end
    fprintf('\n');
    fprintf('  Example (segment #%d at %.1fh):\n', high_beta_idx(1), time_hours(high_beta_idx(1)));
    fprintf('    Delta=%.2f, Theta=%.2f, Alpha=%.2f, Beta=%.2f\n', ...
            delta_values(high_beta_idx(1)), theta_values(high_beta_idx(1)), ...
            alpha_values(high_beta_idx(1)), beta_values(high_beta_idx(1)));
end
fprintf('\n');

%--- Type 3: Paradoxical Pattern (High Delta + High Beta) ---
fprintf('Type 3: Paradoxical Segments (High Delta AND High Beta)\n');
delta_p75 = prctile(delta_values, 75);
high_delta_idx = find(delta_values > delta_p75);
paradoxical_idx = intersect(high_delta_idx, high_beta_idx);

fprintf('  Criteria: Delta > 75th percentile (%.2f) AND Beta > 95th percentile (%.2f)\n', ...
        delta_p75, beta_p95);
fprintf('  Found: %d segments (%.1f%% of total)\n', ...
        length(paradoxical_idx), 100*length(paradoxical_idx)/length(delta_values));

if ~isempty(paradoxical_idx)
    fprintf('  Times: ');
    if length(paradoxical_idx) <= 5
        fprintf('%.1fh ', time_hours(paradoxical_idx));
    else
        fprintf('%.1fh, %.1fh, %.1fh ... (showing first 3)', ...
                time_hours(paradoxical_idx(1)), time_hours(paradoxical_idx(2)), time_hours(paradoxical_idx(3)));
    end
    fprintf('\n');
    fprintf('  Example (segment #%d at %.1fh):\n', paradoxical_idx(1), time_hours(paradoxical_idx(1)));
    fprintf('    Delta=%.2f, Theta=%.2f, Alpha=%.2f, Beta=%.2f\n', ...
            delta_values(paradoxical_idx(1)), theta_values(paradoxical_idx(1)), ...
            alpha_values(paradoxical_idx(1)), beta_values(paradoxical_idx(1)));
    fprintf('  Note: This pattern may indicate:\n');
    fprintf('    - Sleep stage transition\n');
    fprintf('    - Movement artifact\n');
    fprintf('    - REM sleep with delta intrusion\n');
end
fprintf('\n');

%--- Type 4: Extreme Variability (Sudden Changes) ---
fprintf('Type 4: High Variability Segments (Sudden Power Changes)\n');

% Calculate differences between consecutive segments
delta_diff = abs(diff(delta_values));
beta_diff = abs(diff(beta_values));

% Find segments with large changes (>2 SD)
delta_change_threshold = mean(delta_diff) + 2*std(delta_diff);
beta_change_threshold = mean(beta_diff) + 2*std(beta_diff);

high_delta_change_idx = find(delta_diff > delta_change_threshold) + 1; % +1 because diff reduces length
high_beta_change_idx = find(beta_diff > beta_change_threshold) + 1;

fprintf('  Delta sudden changes: %d segments (%.1f%% of total)\n', ...
        length(high_delta_change_idx), 100*length(high_delta_change_idx)/(length(delta_values)-1));
fprintf('  Beta sudden changes: %d segments (%.1f%% of total)\n', ...
        length(high_beta_change_idx), 100*length(high_beta_change_idx)/(length(beta_values)-1));
fprintf('\n');

%% 3. Create Comprehensive Visualizations
fprintf('========================================\n');
fprintf('CREATING VISUALIZATIONS\n');
fprintf('========================================\n\n');

%--- FIGURE 1: Overview of All Unusual Segments ---
figure('Position', [100, 100, 1600, 900], 'Name', 'Unusual Segments Overview');

% Panel 1: Low Delta Segments
subplot(4, 1, 1)
plot(time_hours, delta_values, 'b-', 'LineWidth', 1.5)
hold on
plot(time_hours(low_delta_idx), delta_values(low_delta_idx), 'ro', ...
     'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r')
yline(delta_p05, '--r', 'LineWidth', 2, 'Label', '5th Percentile')
hold off
ylabel('Delta Power (μV²/Hz)', 'FontSize', 11)
title(sprintf('Type 1: Low Delta Segments (%d identified) - Insufficient Deep Sleep', length(low_delta_idx)), ...
      'FontWeight', 'bold', 'FontSize', 12)
legend('Delta Power', 'Low Delta Segments', 'Location', 'best')
grid on
xlim([0 max(time_hours)])

% Panel 2: High Beta Segments
subplot(4, 1, 2)
plot(time_hours, beta_values, 'g-', 'LineWidth', 1.5)
hold on
plot(time_hours(high_beta_idx), beta_values(high_beta_idx), 'ro', ...
     'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r')
yline(beta_p95, '--r', 'LineWidth', 2, 'Label', '95th Percentile')
hold off
ylabel('Beta Power (μV²/Hz)', 'FontSize', 11)
title(sprintf('Type 2: High Beta Segments (%d identified) - Excessive Arousal', length(high_beta_idx)), ...
      'FontWeight', 'bold', 'FontSize', 12)
legend('Beta Power', 'High Beta Segments', 'Location', 'best')
grid on
xlim([0 max(time_hours)])

% Panel 3: Paradoxical Segments
subplot(4, 1, 3)
plot(time_hours, delta_values, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Delta')
hold on
plot(time_hours, beta_values, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Beta')
if ~isempty(paradoxical_idx)
    scatter(time_hours(paradoxical_idx), delta_values(paradoxical_idx), 100, 'r', ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Paradoxical')
end
hold off
ylabel('Power (μV²/Hz)', 'FontSize', 11)
title(sprintf('Type 3: Paradoxical Segments (%d identified) - High Delta + High Beta', length(paradoxical_idx)), ...
      'FontWeight', 'bold', 'FontSize', 12)
legend('Location', 'best')
grid on
xlim([0 max(time_hours)])

% Panel 4: Sudden Changes
subplot(4, 1, 4)
plot(time_hours(2:end), delta_diff, 'b-', 'LineWidth', 1, 'DisplayName', 'Delta Changes')
hold on
plot(time_hours(2:end), beta_diff, 'g-', 'LineWidth', 1, 'DisplayName', 'Beta Changes')
yline(delta_change_threshold, '--b', 'LineWidth', 2, 'DisplayName', 'Delta Threshold')
yline(beta_change_threshold, '--g', 'LineWidth', 2, 'DisplayName', 'Beta Threshold')
plot(time_hours(high_delta_change_idx), delta_diff(high_delta_change_idx-1), 'bo', ...
     'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'b')
plot(time_hours(high_beta_change_idx), beta_diff(high_beta_change_idx-1), 'go', ...
     'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'g')
hold off
xlabel('Time (hours)', 'FontSize', 11)
ylabel('Power Change (μV²/Hz)', 'FontSize', 11)
title(sprintf('Type 4: Sudden Changes (Delta: %d, Beta: %d)', ...
      length(high_delta_change_idx), length(high_beta_change_idx)), ...
      'FontWeight', 'bold', 'FontSize', 12)
legend('Location', 'best')
grid on
xlim([0 max(time_hours)])

sgtitle(sprintf('Unusual Segment Detection - Channel: %s', eeg_labels{primary_channel}), ...
        'FontSize', 16, 'FontWeight', 'bold');

fprintf('✓ Overview figure created\n');



%--- FIGURE 4: Summary Statistics ---
figure('Position', [250, 250, 1200, 800], 'Name', 'Summary Statistics');

% Panel 1: Percentage of unusual segments
subplot(2, 2, 1)
percentages = [100*length(low_delta_idx)/length(delta_values), ...
               100*length(high_beta_idx)/length(beta_values), ...
               100*length(paradoxical_idx)/length(delta_values)];
bar_handle = bar(percentages, 'FaceColor', 'flat');
bar_handle.CData(1,:) = [0.2 0.4 0.8];
bar_handle.CData(2,:) = [0.2 0.8 0.4];
bar_handle.CData(3,:) = [0.8 0.2 0.2];
set(gca, 'XTickLabel', {'Low Delta', 'High Beta', 'Paradoxical'})
ylabel('Percentage of Total Segments', 'FontSize', 11)
title('Unusual Segment Prevalence', 'FontWeight', 'bold')
grid on
ylim([0 max(percentages)*1.2])

% Panel 2: Power distributions for unusual segments
subplot(2, 2, 2)
boxplot([delta_values(low_delta_idx); delta_values], ...
        [ones(length(low_delta_idx),1); 2*ones(length(delta_values),1)], ...
        'Labels', {'Low Delta Segments', 'All Segments'}, 'Colors', 'b')
ylabel('Delta Power (μV²/Hz)', 'FontSize', 11)
title('Low Delta Segment Characteristics', 'FontWeight', 'bold')
grid on

subplot(2, 2, 3)
boxplot([beta_values(high_beta_idx); beta_values], ...
        [ones(length(high_beta_idx),1); 2*ones(length(beta_values),1)], ...
        'Labels', {'High Beta Segments', 'All Segments'}, 'Colors', 'g')
ylabel('Beta Power (μV²/Hz)', 'FontSize', 11)
title('High Beta Segment Characteristics', 'FontWeight', 'bold')
grid on

% Panel 4: All bands for paradoxical segments
if ~isempty(paradoxical_idx)
    subplot(2, 2, 4)
    paradox_powers = [mean(delta_values(paradoxical_idx)), ...
                      mean(theta_values(paradoxical_idx)), ...
                      mean(alpha_values(paradoxical_idx)), ...
                      mean(beta_values(paradoxical_idx))];
    all_powers = [mean(delta_values), mean(theta_values), ...
                  mean(alpha_values), mean(beta_values)];
    
    bar_data = [paradox_powers; all_powers]';
    bar_handle = bar(bar_data);
    bar_handle(1).FaceColor = [0.8 0.2 0.2];
    bar_handle(2).FaceColor = [0.5 0.5 0.5];
    set(gca, 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta'})
    ylabel('Mean Power (μV²/Hz)', 'FontSize', 11)
    title('Paradoxical Segments - All Bands', 'FontWeight', 'bold')
    legend('Paradoxical', 'All Segments', 'Location', 'best')
    grid on
end

sgtitle('Unusual Segment Summary Statistics', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('✓ Summary statistics figure created\n\n');

%% 4. Save Results
fprintf('========================================\n');
fprintf('SAVING RESULTS\n');
fprintf('========================================\n\n');

% Create results structure
results = struct();
results.low_delta.indices = low_delta_idx;
results.low_delta.times = time_hours(low_delta_idx);
results.low_delta.percentage = 100*length(low_delta_idx)/length(delta_values);
results.low_delta.threshold = delta_p05;

results.high_beta.indices = high_beta_idx;
results.high_beta.times = time_hours(high_beta_idx);
results.high_beta.percentage = 100*length(high_beta_idx)/length(beta_values);
results.high_beta.threshold = beta_p95;

results.paradoxical.indices = paradoxical_idx;
results.paradoxical.times = time_hours(paradoxical_idx);
results.paradoxical.percentage = 100*length(paradoxical_idx)/length(delta_values);
results.paradoxical.delta_threshold = delta_p75;
results.paradoxical.beta_threshold = beta_p95;

results.sudden_changes.delta_indices = high_delta_change_idx;
results.sudden_changes.beta_indices = high_beta_change_idx;
results.sudden_changes.delta_threshold = delta_change_threshold;
results.sudden_changes.beta_threshold = beta_change_threshold;

save('../results/step4_unusual_segments.mat', 'results');
fprintf('✓ Results saved to: ../results/step4_unusual_segments.mat\n');

% Save figures
saveas(figure(1), '../figures/step4_unusual_segments_overview.png');
saveas(figure(4), '../figures/step4_summary_statistics.png');
fprintf('✓ All figures saved to figures/ folder\n\n');

%% 5. Final Summary
fprintf('========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n\n');

fprintf('Total Segments Analyzed: %d (%.1f hours)\n', length(time_hours), max(time_hours));
fprintf('Channel: %s\n\n', eeg_labels{primary_channel});

fprintf('UNUSUAL SEGMENTS DETECTED:\n\n');

fprintf('1. Low Delta (Insufficient Deep Sleep):\n');
fprintf('   - Count: %d (%.1f%% of total)\n', length(low_delta_idx), ...
        100*length(low_delta_idx)/length(delta_values));
fprintf('   - Threshold: < %.2f μV²/Hz\n\n', delta_p05);

fprintf('2. High Beta (Excessive Arousal):\n');
fprintf('   - Count: %d (%.1f%% of total)\n', length(high_beta_idx), ...
        100*length(high_beta_idx)/length(beta_values));
fprintf('   - Threshold: > %.2f μV²/Hz\n\n', beta_p95);

fprintf('3. Paradoxical (High Delta + High Beta):\n');
fprintf('   - Count: %d (%.1f%% of total)\n', length(paradoxical_idx), ...
        100*length(paradoxical_idx)/length(delta_values));
fprintf('   - Possible causes: transitions, artifacts, REM intrusions\n\n');

fprintf('4. Sudden Power Changes:\n');
fprintf('   - Delta changes: %d segments\n', length(high_delta_change_idx));
fprintf('   - Beta changes: %d segments\n\n', length(high_beta_change_idx));

fprintf('========================================\n');
fprintf('✓✓✓ STEP 4 COMPLETE! ✓✓✓\n');
fprintf('========================================\n\n');

fprintf('Review the 2 figures created for detailed visual analysis.\n\n');