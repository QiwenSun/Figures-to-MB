%% Fig3 and Fig4 in main text - Combined 3x2 layout with time difference analysis
% Enhanced for high-quality output
% Modified: Add back mean lines, curve size 2.0, remove [88,91] shadow, add (A)-(F) labels

%% ============================================================
% Comparative analysis of gene A and B transcription levels
% v2=0, analyze ω=0.5, 1.0, 1.5, 2.0 for gene A and B in [80, 100] time period
% Enhanced shadow contrast, amplitude-frequency analysis
% Combined 3x2 figure layout with time difference calculation
% High-quality output with anti-aliasing
% ============================================================

%% Clear workspace
clear all; close all; clc;

%% Main program: Compare ODE solutions at different frequencies
fprintf('Starting gene expression system ODE simulation (v2=0)...\n\n');

% Initial conditions [p0, p1, mA, q0, q1, mB]
Y0 = [1; 0; 0; 1; 0; 0];
tspan = [0 100];  % Time range: 0 to 100 hours

% Define frequencies to test
omega_values = [0.5, 1.0, 1.5, 2.0];
n_omega = length(omega_values);

% Pre-store results
segment_data = cell(n_omega, 1);  % Store [80,100] interval data

%% Calculate different frequency cases
for i = 1:n_omega
    omega = omega_values(i);
    fprintf('Calculating ω = %.1f (v2=0)...\n', omega);
    
    % Solve ODE (v2=0)
    [t, Y] = ode45(@(t,Y) myODE_v2_zero(t,Y,omega), tspan, Y0);
    
    % Extract [80, 100] time period data
    idx = find(t >= 80 & t <= 100);
    
    if isempty(idx)
        error('No data found in [80,100] time interval, please check time range settings');
    end
    
    segment_data{i}.time = t(idx);
    segment_data{i}.mA = Y(idx, 3);  % Gene A expression
    segment_data{i}.mB = Y(idx, 6);  % Gene B expression
    
    fprintf('  [80,100] interval data points: %d\n', length(idx));
    fprintf('  Time range: %.2f to %.2f hours\n', min(segment_data{i}.time), max(segment_data{i}.time));
end

fprintf('\nODE solution complete!\n\n');

%% ================ Create high-quality combined 3x2 figure ================
fig = figure('Position', [50, 50, 1800, 1200], ...  % Increased size
             'Name', 'Gene Expression Analysis - Combined 3x2 Layout', ...
             'Color', 'white', ...
             'Renderer', 'painters', ...  % Vector renderer for high quality
             'PaperPositionMode', 'auto', ...
             'InvertHardcopy', 'off', ...  % Keep screen colors
             'PaperUnits', 'inches', ...
             'PaperSize', [12, 8]);        % Paper size for printing

% Define colors - consistent with original
colors = {[0,0,0], [1,0,0], [0,0.5,0], [0,0,1]};  % Black, Red, Green, Blue
shadow_colors_A = {[0.8,0.8,0.8], [1,0.8,0.8], [0.8,1,0.8], [0.8,0.8,1]};  % Light gray, light red, light green, light blue
shadow_colors_B = {[0.6,0.6,0.6], [0.9,0.6,0.6], [0.6,0.9,0.6], [0.6,0.6,0.9]};  % Darker gray, red, green, blue

% Optimized color scheme for amplitude plots
color_alpha = [0, 0.4470, 0.7410];  % MATLAB default blue
color_beta = [0.8500, 0.3250, 0.0980];  % MATLAB default orange
color_ratio = [0.4940, 0.1840, 0.5560];  % MATLAB default purple

%% ================ Subplot 1-4: Time series for four frequencies ================
for i = 1:4
    subplot(3, 2, i);
    
    % First draw shadow areas (ensure below curves)
    % Calculate statistics
    mA_mean = mean(segment_data{i}.mA);
    mA_std = std(segment_data{i}.mA);
    mA_min = min(segment_data{i}.mA);
    mA_max = max(segment_data{i}.mA);
    mA_amplitude = (mA_max - mA_min)/2;
    
    mB_mean = mean(segment_data{i}.mB);
    mB_std = std(segment_data{i}.mB);
    mB_min = min(segment_data{i}.mB);
    mB_max = max(segment_data{i}.mB);
    mB_amplitude = (mB_max - mB_min)/2;
    
    % ==== Calculate time difference in specific interval [88, 91] ====
    time_data = segment_data{i}.time;
    mA_data = segment_data{i}.mA;
    mB_data = segment_data{i}.mB;
    
    % Find indices within [88, 91] time interval
    idx_interval = find(time_data >= 88 & time_data <= 91);
    
    if ~isempty(idx_interval)
        % Extract data in the specific interval
        time_interval = time_data(idx_interval);
        mA_interval = mA_data(idx_interval);
        mB_interval = mB_data(idx_interval);
        
        % Find time of maximum for gene A in this interval
        [mA_max_val, idx_mA_max] = max(mA_interval);
        t_max_A = time_interval(idx_mA_max);
        
        % Find time of maximum for gene B in this interval
        [mB_max_val, idx_mB_max] = max(mB_interval);
        t_max_B = time_interval(idx_mB_max);
        
        % Calculate time difference (B - A)
        time_diff = t_max_B - t_max_A;
        
        % Add markers to the plot for visualization
        hold on;
        plot(t_max_A, mA_max_val, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
        plot(t_max_B, mB_max_val, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
        
        % Display time difference in the text box
        time_diff_text = sprintf('Δt = %.3f hours', time_diff);
    else
        time_diff_text = 'No data in [88,91] interval';
        time_diff = NaN;
    end
    
    % Store the time difference for later use if needed
    segment_data{i}.time_diff = time_diff;
    segment_data{i}.t_max_A = t_max_A;
    segment_data{i}.t_max_B = t_max_B;
    
    % Draw gene A shadow area
    fill_x = [segment_data{i}.time; flipud(segment_data{i}.time)];
    fill_y_A = [ones(size(segment_data{i}.mA))*(mA_mean - mA_std); 
                flipud(ones(size(segment_data{i}.mA))*(mA_mean + mA_std))];
    h_fill_A = fill(fill_x, fill_y_A, shadow_colors_A{i}, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    
    hold on;
    
    % Draw gene B shadow area
    fill_y_B = [ones(size(segment_data{i}.mB))*(mB_mean - mB_std); 
                flipud(ones(size(segment_data{i}.mB))*(mB_mean + mB_std))];
    h_fill_B = fill(fill_x, fill_y_B, shadow_colors_B{i}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Then draw curves (above shadows) - capture handles for legend
    h1 = plot(segment_data{i}.time, segment_data{i}.mA, ...
         'Color', colors{i}, 'LineWidth', 2.0);  % MODIFIED: LineWidth = 2.0
    
    h2 = plot(segment_data{i}.time, segment_data{i}.mB, ...
         'Color', colors{i}, 'LineWidth', 2.0, 'LineStyle', '--');  % MODIFIED: LineWidth = 2.0
    
    % ADDED BACK: Mean lines with thinner lines
    h_mean_A = plot([80, 100], [mA_mean, mA_mean], 'k-', 'LineWidth', 1.2, 'LineStyle', ':');
    h_mean_B = plot([80, 100], [mB_mean, mB_mean], 'k-', 'LineWidth', 1.2, 'LineStyle', '-.');
    
    % Set subplot properties
    xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Transcription Level', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add panel labels with (A)(B)(C)(D) prefix
    if i == 1
        title('(A) \omega = 0.5', 'FontSize', 14, 'FontWeight', 'bold');
    elseif i == 2
        title('(B) \omega = 1.0', 'FontSize', 14, 'FontWeight', 'bold');
    elseif i == 3
        title('(C) \omega = 1.5', 'FontSize', 14, 'FontWeight', 'bold');
    elseif i == 4
        title('(D) \omega = 2.0', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Add statistical information text (repositioned for better visibility)
    % Position the text box at top right
    y_max = max([max(segment_data{i}.mA), max(segment_data{i}.mB)]);
    
    % Create text box with statistical info including time difference
    text_content = sprintf(['Gene α:\nMean: %.3f\nAmp: %.3f\n' ...
                           'Gene β:\nMean: %.3f\nAmp: %.3f\n' ...
                           '%s'], ...
                 mA_mean, mA_amplitude, mB_mean, mB_amplitude, time_diff_text);
    
    % Position at top right corner
    h_text = text(98, y_max*0.92+1.5, text_content, ...
         'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    
    % Add legend only in first subplot
    if i == 1
        h_legend = legend([h1, h2], {'m_α(t)', 'm_β(t)'}, ...
               'Location', 'best', 'FontSize', 10);
        set(h_legend, 'Box', 'off');
    end
    
    grid on;
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on', ...
             'GridLineStyle', ':', 'GridAlpha', 0.3);
    xlim([80, 100]);
    
    % Set y-axis limits to be consistent across panels
    y_limits = [min([mA_min, mB_min])*0.95, max([mA_max, mB_max])*1.05];
    ylim(y_limits);
    
    % REMOVED: [88, 91] shaded background - commented out
    % yl = ylim;
    % h_highlight = fill([88, 91, 91, 88], [yl(1), yl(1), yl(2), yl(2)], ...
    %      [0.9, 0.9, 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Simplified layer management - control by drawing order
    % Draw shadows first, then curves, finally text
    
    % Set drawing order - ensure text is on top
    % Use uistack but in a safer way
    try
        uistack(h_text, 'top');
    catch
        % If uistack fails, ignore error
    end
end

%% ================ Subplot 5: Amplitude-Frequency Relationship ================
subplot(3, 2, 5);
hold on;

% Prepare statistical data for amplitude plots
stat_amplitudes_A = zeros(n_omega, 1);
stat_amplitudes_B = zeros(n_omega, 1);
time_diffs = zeros(n_omega, 1);  % Store time differences for frequency plot

for i = 1:n_omega
    % Gene A statistics
    mA_data = segment_data{i}.mA;
    stat_amplitudes_A(i) = (max(mA_data) - min(mA_data))/2;
    
    % Gene B statistics
    mB_data = segment_data{i}.mB;
    stat_amplitudes_B(i) = (max(mB_data) - min(mB_data))/2;
    
    % Time difference
    time_diffs(i) = segment_data{i}.time_diff;
end

% Plot gene α amplitude-frequency relationship
h_alpha = plot(omega_values, stat_amplitudes_A, 'o-', ...
     'Color', color_alpha, ...
     'LineWidth', 2.0, 'MarkerSize', 10, 'MarkerFaceColor', color_alpha);  % MODIFIED: LineWidth = 2.0

% Plot gene β amplitude-frequency relationship
h_beta = plot(omega_values, stat_amplitudes_B, 's--', ...
     'Color', color_beta, ...
     'LineWidth', 2.0, 'MarkerSize', 10, 'MarkerFaceColor', color_beta);  % MODIFIED: LineWidth = 2.0

xlabel('Frequency \omega (rad/hour)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
title('(E) Amplitude-Frequency Relationship', 'FontSize', 14, 'FontWeight', 'bold');  % MODIFIED: Added (E)

h_legend_e = legend([h_alpha, h_beta], {'Gene \alpha', 'Gene \beta'}, ...
               'Location', 'best', 'FontSize', 11);
set(h_legend_e, 'Box', 'off');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on', ...
         'GridLineStyle', ':', 'GridAlpha', 0.3);

% Add numerical labels
for i = 1:n_omega
    % Position labels to avoid overlap
    if i == 1
        offset_x = 0.15;
        offset_y_A = 0.015;
        offset_y_B = -0.015;
    elseif i == 2
        offset_x = 0.15;
        offset_y_A = 0.015;
        offset_y_B = -0.015;
    elseif i == 3
        offset_x = 0.15;
        offset_y_A = -0.02;
        offset_y_B = 0.02;
    else
        offset_x = 0.15;
        offset_y_A = -0.02;
        offset_y_B = 0.02;
    end
    
    % Add text labels
    text(omega_values(i)+offset_x, stat_amplitudes_A(i)+offset_y_A, ...
         sprintf('%.3f', stat_amplitudes_A(i)), ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', color_alpha);
    
    text(omega_values(i)+offset_x, stat_amplitudes_B(i)+offset_y_B, ...
         sprintf('%.3f', stat_amplitudes_B(i)), ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', color_beta);
end

%% ================ Subplot 6: Amplitude Ratio (β/α) and Time Difference ================
subplot(3, 2, 6);
hold on;

% Calculate amplitude ratio
amplitude_ratio = stat_amplitudes_B ./ stat_amplitudes_A;

% Plot amplitude ratio on left y-axis
yyaxis left;
h_ratio = plot(omega_values, amplitude_ratio, 'd-', ...
     'Color', color_ratio, ...
     'LineWidth', 2.0, 'MarkerSize', 10, 'MarkerFaceColor', color_ratio);  % MODIFIED: LineWidth = 2.0

ylabel('Amplitude Ratio', 'FontSize', 12, 'FontWeight', 'bold', 'Color', color_ratio);
ylim([0, 1.2]);  % Adjust based on data
set(gca, 'YColor', color_ratio);

% Plot time difference on right y-axis
yyaxis right;
h_time_diff = plot(omega_values, time_diffs, '^-', ...
     'Color', [0.2, 0.6, 0.2], ...
     'LineWidth', 2.0, 'MarkerSize', 10, 'MarkerFaceColor', [0.2, 0.6, 0.2]);  % MODIFIED: LineWidth = 2.0

ylabel('Time Difference', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2, 0.6, 0.2]);
set(gca, 'YColor', [0.2, 0.6, 0.2]);

xlabel('Frequency \omega (rad/hour)', 'FontSize', 12, 'FontWeight', 'bold');
title('(F) Amplitude Ratio & Time Difference', 'FontSize', 14, 'FontWeight', 'bold');  % MODIFIED: Added (F)

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on', ...
         'GridLineStyle', ':', 'GridAlpha', 0.3);

% REMOVED: Reference line at y=1 for amplitude ratio
% yyaxis left;
% h_ref_ratio = plot(xlim, [1, 1], 'k--', 'LineWidth', 2.0, 'Color', [0.5, 0.5, 0.5]);

% RETAINED: Reference line at y=0 for time difference
yyaxis right;
h_ref_time = plot(xlim, [0, 0], '--', 'LineWidth', 2.0, 'Color', [0.5, 0.6, 0.5]);  % MODIFIED: LineWidth = 2.0

% Add legend - REMOVED: 'Ratio = 1' from legend
h_legend_f = legend([h_ratio, h_time_diff], ...
       {'Amplitude ratio', 'Time diff (\beta-\alpha)'}, ...
       'Location', 'southwest', 'FontSize', 11);
set(h_legend_f, 'Box', 'off');

% Add amplitude ratio labels (left y-axis)
yyaxis left;
for i = 1:n_omega
    text(omega_values(i)+0.05, amplitude_ratio(i)-0.1, ...
         sprintf('%.3f', amplitude_ratio(i)), ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', color_ratio);
end

% Add time difference labels (right y-axis)
yyaxis right;
for i = 1:n_omega
    text(omega_values(i)+0.05, time_diffs(i)+0.1, ...
         sprintf('%.3f', time_diffs(i)), ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', [0.2, 0.6, 0.2]);
end

%% ================ Save High-Quality Figures ================
fprintf('\nSaving high-quality figures...\n');

% Ensure figure is up to date
drawnow;

% Save in multiple formats for different uses
try
    % 1. PNG format for general use (high resolution)
    print(fig, '-dpng', '-r600', 'gene_expression_combined_hires.png');
    fprintf('  High-res PNG saved: gene_expression_combined_hires.png (600 DPI)\n');
    
    % 2. PDF format for publications (vector format, infinite zoom)
    print(fig, '-dpdf', '-r600', 'gene_expression_combined.pdf');
    fprintf('  PDF saved: gene_expression_combined.pdf (vector format)\n');
    
    % 3. SVG format for editing (vector format)
    print(fig, '-dsvg', '-r600', 'gene_expression_combined.svg');
    fprintf('  SVG saved: gene_expression_combined.svg (editable vector)\n');
    
    % 4. EPS format for LaTeX
    print(fig, '-depsc2', '-r600', 'gene_expression_combined.eps');
    fprintf('  EPS saved: gene_expression_combined.eps (for LaTeX)\n');
    
    % 5. TIFF format for journals requiring it
    print(fig, '-dtiff', '-r600', 'gene_expression_combined.tiff');
    fprintf('  TIFF saved: gene_expression_combined.tiff (600 DPI)\n');
    
    % 6. FIG format for MATLAB editing
    savefig(fig, 'gene_expression_combined.fig');
    fprintf('  FIG saved: gene_expression_combined.fig (MATLAB format)\n');
    
    fprintf('\nAll figures saved successfully!\n');
    fprintf('Recommendations:\n');
    fprintf('  - Use PDF/SVG for publications (vector, no blur when zooming)\n');
    fprintf('  - Use PNG for presentations/web\n');
    fprintf('  - Use FIG to edit in MATLAB later\n');
    
catch ME
    fprintf('Error saving figures: %s\n', ME.message);
    fprintf('Trying alternative save methods...\n');
    
    % Try alternative save methods
    try
        saveas(fig, 'gene_expression_combined.png');
        fprintf('  Saved using saveas: gene_expression_combined.png\n');
    catch
        fprintf('  Could not save figures. Please check write permissions.\n');
    end
end

%% ================ Display Summary Statistics ================
fprintf('\n=== Summary Statistics ===\n');
fprintf('Frequency | Gene α Amplitude | Gene β Amplitude | Ratio (β/α) | Time Diff (hours)\n');
fprintf('----------|------------------|------------------|-------------|-------------------\n');
for i = 1:n_omega
    fprintf('   %.1f    |      %.4f      |      %.4f      |    %.4f    |      %.4f\n', ...
            omega_values(i), stat_amplitudes_A(i), stat_amplitudes_B(i), ...
            amplitude_ratio(i), time_diffs(i));
end

%% ================ ODE Function Definition (v2=0) ================
function dYdt = myODE_v2_zero(t, Y, omega)
    % Unified gene expression system ODE function (v2=0)
    % Input parameters:
    %   t - time
    %   Y - state variables [p0; p1; mA; q0; q1; mB]
    %   omega - driving frequency parameter
    
    % Unpack state variables
    p0 = Y(1);  % Gene A promoter state 0
    p1 = Y(2);  % Gene A promoter state 1
    mA = Y(3);  % Gene A mRNA level
    q0 = Y(4);  % Gene B promoter state 0
    q1 = Y(5);  % Gene B promoter state 1
    mB = Y(6);  % Gene B mRNA level

    % System parameters - Gene A
    k_on_A = 0.5;    % Gene A ON rate
    k_off_A = 0.5;   % Gene A OFF rate
    gamma_A = 1;     % Gene A mRNA degradation rate
    
    % Time-varying transcription rate - includes cosine oscillation
    a = 5;      % oscillation amplitude of Gene A
    V0_A = 20;  % Basal synthesis rate of Gene A
    v1 = V0_A + a * cos(omega * t);
    
    % System parameters - Gene B (v2=0!)
    k_on_B = 0.5;    % Gene B ON rate
    k_off_B = 0.5;   % Gene B OFF rate
    gamma_B = 1;     % Gene B mRNA degradation rate 
    V0_B = 0;        % Gene B Basal transcription rate = 0 !!!
    k_B = 2;         % Gene A regulation strength on Gene B
    
    % Define differential equations
    % Gene A promoter dynamics
    dp0dt = k_on_A * p1 - k_off_A * p0;      % p0 rate of change
    dp1dt = k_off_A * p0 - k_on_A * p1;      % p1 rate of change
    
    % Gene A mRNA dynamics
    dmAdt = -gamma_A * mA + v1 * p1;     % mA rate of change
    
    % Gene B promoter dynamics
    dq0dt = k_on_B * q1 - k_off_B * q0;      % q0 rate of change
    dq1dt = k_off_B * q0 - k_on_B * q1;      % q1 rate of change
    
    % Gene B mRNA dynamics (regulated by Gene A, and v2=0)
    dmBdt = -gamma_B * mB + (V0_B + k_B * mA) * q1;  % Key modification: v2=0
    
    % Package derivatives into vector
    dYdt = [dp0dt; dp1dt; dmAdt; dq0dt; dq1dt; dmBdt];
end