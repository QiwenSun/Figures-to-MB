%% Fig7 in main text - Modified for academic journal publication
function gillespie_gene_expression()
    % Gillespie Stochastic Simulation Algorithm for gene A and gene B transcription
    % Genes randomly switch between ON and OFF states
    % mRNA synthesis only occurs in ON state
    % Modified for academic journal publication standards
    
    clear; clc; close all;
    
    %% Define academic publication color scheme
    % Use Nature/Science style colors (distinguishable in grayscale)
    color_geneA = [0, 0.45, 0.74];     % Blue - distinguishable
    color_geneB = [0.85, 0.33, 0.10];  % Orange/Red - good contrast
    color_geneA_light = [0.7, 0.85, 1.0]; % Light blue for fills
    color_geneB_light = [1.0, 0.8, 0.7];  % Light orange for fills
    
    % Line styles for better grayscale distinction
    line_style_geneA = '-';  % Solid for Gene A
    line_style_geneB = '--'; % Dashed for Gene B
    
    % Marker styles
    marker_max = 'o';      % Circle for maxima
    marker_min = 's';      % Square for minima
    
    % Set publication font properties
    font_name = 'Helvetica';  % Common publication font
    font_size_title = 12;
    font_size_axis = 11;
    font_size_legend = 10;
    font_size_text = 9;
    
    %% Parameter settings
    num_simulations = 10000;  % Number of simulation repetitions
    T_max = 100;              % Total simulation time
    dt = 0.1;                 % Time step (for recording results)
    
    % Gene state switching parameters
    k_on_A = 0.5;         % Gene A ON rate
    k_off_A = 0.5;        % Gene A OFF rate
    k_on_B = 0.5;         % Gene B ON rate
    k_off_B = 0.5;        % Gene B OFF rate
    
    % mRNA degradation parameters
    gamma_A = 1;          % Gene A mRNA degradation rate
    gamma_B = 1;          % Gene B mRNA degradation rate
    
    % Gene A synthesis parameters (periodic)
    V0_A = 20;            % Basal synthesis rate
    a = 10;               % Amplitude
    omega = 1;            % Angular frequency
    
    % Gene B synthesis parameters (proportional to gene A mRNA amount)
    k_B = 2;              % Proportionality coefficient
    
    %% Initialize storage arrays
    time_points = 0:dt:T_max;
    num_time_points = length(time_points);
    
    % Store results for each simulation
    mRNA_A_all = zeros(num_simulations, num_time_points);
    mRNA_B_all = zeros(num_simulations, num_time_points);
    state_A_all = zeros(num_simulations, num_time_points);
    state_B_all = zeros(num_simulations, num_time_points);
    
    %% Perform multiple simulations
    fprintf('Running Gillespie simulation for publication figure...\n');
    
    for sim = 1:num_simulations
        if mod(sim, 1000) == 0
            fprintf('Completed %d/%d simulations\n', sim, num_simulations);
        end
        
        % Initialize variables for single simulation
        t_current = 0;
        mRNA_A = 0;
        mRNA_B = 0;
        state_A = 0;  % 0=OFF, 1=ON
        state_B = 0;  % 0=OFF, 1=ON
        
        time_index = 1;
        mRNA_A_all(sim, time_index) = mRNA_A;
        mRNA_B_all(sim, time_index) = mRNA_B;
        state_A_all(sim, time_index) = state_A;
        state_B_all(sim, time_index) = state_B;
        
        % Gillespie algorithm main loop
        while t_current < T_max
            % Calculate all possible reaction rates
            rates = zeros(6, 1);
            
            % 1: Gene A state switching
            if state_A == 0
                rates(1) = k_on_A;  % OFF -> ON
            else
                rates(1) = k_off_A; % ON -> OFF
            end
            
            % 2: Gene B state switching
            if state_B == 0
                rates(2) = k_on_B;  % OFF -> ON
            else
                rates(2) = k_off_B; % ON -> OFF
            end
            
            % 3: Gene A mRNA synthesis (only in ON state)
            if state_A == 1
                V_A = V0_A + a * cos(omega * t_current);
                rates(3) = V_A;
            else
                rates(3) = 0;
            end
            
            % 4: Gene B mRNA synthesis (only in ON state and proportional to gene A mRNA)
            if state_B == 1
                V_B = k_B * mRNA_A;
                rates(4) = V_B;
            else
                rates(4) = 0;
            end
            
            % 5: Gene A mRNA degradation
            rates(5) = gamma_A * mRNA_A;
            
            % 6: Gene B mRNA degradation
            rates(6) = gamma_B * mRNA_B;
            
            % Calculate total reaction rate
            total_rate = sum(rates);
            
            if total_rate == 0
                % If no reactions are possible, break the loop
                break;
            end
            
            % Generate two random numbers
            r1 = rand();
            r2 = rand();
            
            % Calculate time to next reaction
            tau = -log(r1) / total_rate;
            
            % Determine which reaction occurs
            cumulative_rate = 0;
            reaction_index = 0;
            for i = 1:6
                cumulative_rate = cumulative_rate + rates(i);
                if r2 * total_rate <= cumulative_rate
                    reaction_index = i;
                    break;
                end
            end
            
            % Update time
            t_current = t_current + tau;
            
            % Execute reaction
            switch reaction_index
                case 1  % Gene A state switching
                    state_A = 1 - state_A;  % Toggle state
                    
                case 2  % Gene B state switching
                    state_B = 1 - state_B;  % Toggle state
                    
                case 3  % Gene A mRNA synthesis
                    mRNA_A = mRNA_A + 1;
                    
                case 4  % Gene B mRNA synthesis
                    mRNA_B = mRNA_B + 1;
                    
                case 5  % Gene A mRNA degradation
                    mRNA_A = mRNA_A - 1;
                    if mRNA_A < 0
                        mRNA_A = 0;
                    end
                    
                case 6  % Gene B mRNA degradation
                    mRNA_B = mRNA_B - 1;
                    if mRNA_B < 0
                        mRNA_B = 0;
                    end
            end
            
            % Record results
            while t_current >= time_points(time_index) && time_index < num_time_points
                time_index = time_index + 1;
                mRNA_A_all(sim, time_index) = mRNA_A;
                mRNA_B_all(sim, time_index) = mRNA_B;
                state_A_all(sim, time_index) = state_A;
                state_B_all(sim, time_index) = state_B;
            end
            
            % Break if exceeding maximum time
            if t_current > T_max
                break;
            end
        end
        
        % Fill remaining time points (if simulation ends early)
        while time_index < num_time_points
            time_index = time_index + 1;
            mRNA_A_all(sim, time_index) = mRNA_A;
            mRNA_B_all(sim, time_index) = mRNA_B;
            state_A_all(sim, time_index) = state_A;
            state_B_all(sim, time_index) = state_B;
        end
    end
    
    %% Calculate averages and statistics
    mRNA_A_mean = mean(mRNA_A_all, 1);
    mRNA_B_mean = mean(mRNA_B_all, 1);
    state_A_mean = mean(state_A_all, 1);
    state_B_mean = mean(state_B_all, 1);
    
    % Calculate standard deviations for error bands
    mRNA_A_std = std(mRNA_A_all, 0, 1);
    mRNA_B_std = std(mRNA_B_all, 0, 1);
    
    %% Create publication-quality figure with optimized layout
    fig = figure('Position', [100, 100, 1800, 1000], ...
                 'Color', 'white', ...
                 'Renderer', 'painters', ...  % Vector renderer for high quality
                 'PaperPositionMode', 'auto', ...
                 'InvertHardcopy', 'off', ...  % Keep screen colors
                 'PaperUnits', 'inches', ...
                 'PaperSize', [15, 10]);       % Larger paper size for better resolution
    
    % Set default font for all axes
    set(groot, 'DefaultAxesFontName', font_name);
    set(groot, 'DefaultTextFontName', font_name);
    
    %% Subplot 1: Single stochastic realization (A)
    subplot(2, 4, 1);
    hold on;
    
    % Plot with distinct line styles - increase line width for EPS output
    h1 = plot(time_points, mRNA_A_all(1, :), line_style_geneA, ...
              'Color', color_geneA, 'LineWidth', 2.0);  % Increased from 1.8
    h2 = plot(time_points, mRNA_B_all(1, :), line_style_geneB, ...
              'Color', color_geneB, 'LineWidth', 2.0);  % Increased from 1.8
    
    % Formatting
    xlabel('Time (hours)', 'FontSize', font_size_axis, 'FontWeight', 'normal');
    ylabel('mRNA copy number', 'FontSize', font_size_axis, 'FontWeight', 'normal');
    title('(A) Single stochastic realization', ...
          'FontSize', font_size_title, 'FontWeight', 'bold');
    
    % Add legend with publication style
    legend([h1, h2], {'Gene α', 'Gene β'}, ...
           'Location', 'best', 'FontSize', font_size_legend, ...
           'Box', 'off', 'EdgeColor', 'none');
    
    % Grid and axis properties
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, ...
             'Box', 'on', 'LineWidth', 1.2, ...  % Increased axis line width
             'FontSize', font_size_axis-1);
    xlim([0, T_max]);
    ylim([0, max([max(mRNA_A_all(1, :)), max(mRNA_B_all(1, :))]) * 1.1]);
    
    %% Subplot 2: Ensemble averages with statistical analysis (B)
    subplot(2, 4, [2, 3]); % Span two columns
    hold on;
    
    % Define analysis time window [67, 80]
    time_mask = time_points >= 67 & time_points <= 80;
    time_analysis = time_points(time_mask);
    mA_mean_analysis = mRNA_A_mean(time_mask);
    mB_mean_analysis = mRNA_B_mean(time_mask);
    mA_std_analysis = mRNA_A_std(time_mask);
    mB_std_analysis = mRNA_B_std(time_mask);
    
    % Plot means with error bands (standard deviation)
    fill_x = [time_analysis, fliplr(time_analysis)];
    
    % Gene A error band - use simpler fill for EPS compatibility
    fill_y_A = [mA_mean_analysis - mA_std_analysis, ...
                fliplr(mA_mean_analysis + mA_std_analysis)];
    h_fill_A = fill(fill_x, fill_y_A, color_geneA, ...
                    'FaceAlpha', 0.25, ...  % Reduced alpha for better EPS output
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Gene B error band
    fill_y_B = [mB_mean_analysis - mB_std_analysis, ...
                fliplr(mB_mean_analysis + mB_std_analysis)];
    h_fill_B = fill(fill_x, fill_y_B, color_geneB, ...
                    'FaceAlpha', 0.25, ...  % Reduced alpha for better EPS output
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Plot mean trajectories with thicker lines for EPS
    h3 = plot(time_analysis, mA_mean_analysis, line_style_geneA, ...
              'Color', color_geneA, 'LineWidth', 3.0);  % Increased from 2.5
    h4 = plot(time_analysis, mB_mean_analysis, line_style_geneB, ...
              'Color', color_geneB, 'LineWidth', 3.0);  % Increased from 2.5
    
    % Find and mark extremal points for one complete cycle
    [max_vals_A, max_idx_A] = findpeaks(mA_mean_analysis);
    [min_vals_A, min_idx_A] = findpeaks(-mA_mean_analysis);
    min_vals_A = -min_vals_A;
    
    [max_vals_B, max_idx_B] = findpeaks(mB_mean_analysis);
    [min_vals_B, min_idx_B] = findpeaks(-mB_mean_analysis);
    min_vals_B = -min_vals_B;
    
    % Select middle cycle for clarity
    if length(max_idx_A) >= 2 && length(min_idx_A) >= 2
        mid_idx = floor(length(max_idx_A)/2);
        time_max_A = time_analysis(max_idx_A(mid_idx));
        val_max_A = max_vals_A(mid_idx);
        
        min_after_max = min_idx_A(min_idx_A > max_idx_A(mid_idx));
        if ~isempty(min_after_max)
            time_min_A = time_analysis(min_after_max(1));
            val_min_A = min_vals_A(min_idx_A == min_after_max(1));
        else
            time_min_A = NaN;
            val_min_A = NaN;
        end
    else
        [val_max_A, idx_max_A] = max(mA_mean_analysis);
        [val_min_A, idx_min_A] = min(mA_mean_analysis);
        time_max_A = time_analysis(idx_max_A);
        time_min_A = time_analysis(idx_min_A);
    end
    
    if length(max_idx_B) >= 2 && length(min_idx_B) >= 2
        mid_idx = floor(length(max_idx_B)/2);
        time_max_B = time_analysis(max_idx_B(mid_idx));
        val_max_B = max_vals_B(mid_idx);
        
        min_after_max = min_idx_B(min_idx_B > max_idx_B(mid_idx));
        if ~isempty(min_after_max)
            time_min_B = time_analysis(min_after_max(1));
            val_min_B = min_vals_B(min_idx_B == min_after_max(1));
        else
            time_min_B = NaN;
            val_min_B = NaN;
        end
    else
        [val_max_B, idx_max_B] = max(mB_mean_analysis);
        [val_min_B, idx_min_B] = min(mB_mean_analysis);
        time_max_B = time_analysis(idx_max_B);
        time_min_B = time_analysis(idx_min_B);
    end
    
    % Mark extremal points with publication-style markers
    % Increase marker size and line width for better EPS output
    plot(time_max_A, val_max_A, marker_max, ...
         'MarkerSize', 12, 'MarkerFaceColor', color_geneA, ...  % Increased from 10
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % Increased from 1.2
    plot(time_min_A, val_min_A, marker_min, ...
         'MarkerSize', 12, 'MarkerFaceColor', color_geneA, ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(time_max_B, val_max_B, marker_max, ...
         'MarkerSize', 12, 'MarkerFaceColor', color_geneB, ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(time_min_B, val_min_B, marker_min, ...
         'MarkerSize', 12, 'MarkerFaceColor', color_geneB, ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Add time labels - increase font size for better EPS output
    text(time_max_A, val_max_A*1.05, sprintf('t=%.1f', time_max_A), ...
         'FontSize', font_size_text+1, 'HorizontalAlignment', 'center', ...  % Increased from font_size_text
         'BackgroundColor', 'white', 'Margin', 1);
    text(time_min_A, val_min_A*0.95, sprintf('t=%.1f', time_min_A), ...
         'FontSize', font_size_text+1, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'Margin', 1);
    text(time_max_B, val_max_B*1.05, sprintf('t=%.1f', time_max_B), ...
         'FontSize', font_size_text+1, 'HorizontalAlignment', 'center', ...
         'BackgroundColor', 'white', 'Margin', 1);
    text(time_min_B, val_min_B*0.95, sprintf('t=%.1f', time_min_B), ...
         'FontSize', font_size_text+1, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'Margin', 1);
    
    % Formatting
    xlabel('Time (hours)', 'FontSize', font_size_axis, 'FontWeight', 'normal');
    ylabel('⟨mRNA⟩ (copy number)', 'FontSize', font_size_axis, 'FontWeight', 'normal');
    title('(B) Ensemble averages (n=10,000)', ...
          'FontSize', font_size_title, 'FontWeight', 'bold');
    
    % Add legend for error bands
    h_fill_A_dummy = fill(NaN, NaN, color_geneA, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h_fill_B_dummy = fill(NaN, NaN, color_geneB, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    
    legend([h3, h_fill_A_dummy, h4, h_fill_B_dummy], ...
           {'⟨m_α⟩', '±σ_α', '⟨m_β⟩', '±σ_β'}, ...
           'Location', 'best', 'FontSize', font_size_legend, ...
           'Box', 'off', 'NumColumns', 2);
    
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, ...
             'Box', 'on', 'LineWidth', 1.2, ...
             'FontSize', font_size_axis-1);
    xlim([67, 80]);
    ylim([min([mA_mean_analysis, mB_mean_analysis])*0.9, ...
          max([mA_mean_analysis, mB_mean_analysis])*1.1]);
    
    %% Subplot 3: Global mRNA distribution (C)
    subplot(2, 4, 4);
    hold on;
    
    % Flatten all simulation data
    mRNA_A_global = mRNA_A_all(:);
    mRNA_B_global = mRNA_B_all(:);
    
    % Calculate statistics
    mean_A_global = mean(mRNA_A_global);
    mean_B_global = mean(mRNA_B_global);
    var_A_global = var(mRNA_A_global);
    var_B_global = var(mRNA_B_global);
    noise_A_global = var_A_global / (mean_A_global^2);
    noise_B_global = var_B_global / (mean_B_global^2);
    
    % Create histogram with publication style
    bin_edges = 0:2:max([mRNA_A_global; mRNA_B_global])+2;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    [counts_A, ~] = histcounts(mRNA_A_global, bin_edges, 'Normalization', 'pdf');
    [counts_B, ~] = histcounts(mRNA_B_global, bin_edges, 'Normalization', 'pdf');
    
    % Plot as bar with increased edge width for EPS
    bar_width = 0.8;
    bar(bin_centers - bar_width/2, counts_A, bar_width, ...
        'FaceColor', color_geneA, 'FaceAlpha', 0.8, ...  % Increased alpha
        'EdgeColor', color_geneA*0.5, 'LineWidth', 1.2);  % Increased line width
    
    bar(bin_centers + bar_width/2, counts_B, bar_width, ...
        'FaceColor', color_geneB, 'FaceAlpha', 0.8, ...  % Increased alpha
        'EdgeColor', color_geneB*0.5, 'LineWidth', 1.2);  % Increased line width
    
    % Add mean lines with increased width
    y_max = max([counts_A, counts_B]) * 1.1;
    line([mean_A_global, mean_A_global], [0, y_max*0.9], ...
         'Color', color_geneA, 'LineWidth', 2.8, 'LineStyle', '--');  % Increased from 2.5
    line([mean_B_global, mean_B_global], [0, y_max*0.9], ...
         'Color', color_geneB, 'LineWidth', 2.8, 'LineStyle', '--');  % Increased from 2.5
    
    % Statistical annotation
    stat_text = {
        sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Gene α:', ...
                color_geneA(1), color_geneA(2), color_geneA(3))
        sprintf('⟨m⟩ = %.2f', mean_A_global)
        sprintf('σ²/⟨m⟩² = %.4f', noise_A_global)
        ''
        sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Gene β:', ...
                color_geneB(1), color_geneB(2), color_geneB(3))
        sprintf('⟨m⟩ = %.2f', mean_B_global)
        sprintf('σ²/⟨m⟩² = %.4f', noise_B_global)
    };
    
    text(0.25, 0.95, stat_text, ...
         'Units', 'normalized', 'FontSize', font_size_text+1, ...  % Increased font size
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'BackgroundColor', 'white', 'EdgeColor', [0.5, 0.5, 0.5], ...
         'Margin', 2, 'Interpreter', 'tex');
    
    % Formatting
    xlabel('mRNA copy number', 'FontSize', font_size_axis, 'FontWeight', 'normal');
    ylabel('Probability density', 'FontSize', font_size_axis, 'FontWeight', 'normal');
    title('(C) Stationary distribution', ...
          'FontSize', font_size_title, 'FontWeight', 'bold');
    
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, ...
             'Box', 'on', 'LineWidth', 1.2, ...
             'FontSize', font_size_axis-1);
    xlim([bin_edges(1)-1, bin_edges(end)]);
    
    %% Subplots 5-8: Time-specific distributions at extremal points (D-G)
    extremal_times = [time_max_A, time_min_A, time_max_B, time_min_B];
    extremal_labels = {'D', 'E', 'F', 'G'};
    extremal_descriptions = {'α maximum', 'α minimum', 'β maximum', 'β minimum'};
    
    for i = 1:4
        subplot(2, 4, 4+i);
        hold on;
        
        % Get data at specific time point
        time_val = extremal_times(i);
        time_idx = find(abs(time_points - time_val) < dt/2, 1);
        
        if ~isempty(time_idx)
            mRNA_A_time = mRNA_A_all(:, time_idx);
            mRNA_B_time = mRNA_B_all(:, time_idx);
            
            % Calculate statistics
            mean_A_time = mean(mRNA_A_time);
            mean_B_time = mean(mRNA_B_time);
            var_A_time = var(mRNA_A_time);
            var_B_time = var(mRNA_B_time);
            noise_A_time = var_A_time / (mean_A_time^2);
            noise_B_time = var_B_time / (mean_B_time^2);
            
            % Create histogram
            min_val = min([mRNA_A_time; mRNA_B_time]);
            max_val = max([mRNA_A_time; mRNA_B_time]);
            bin_edges_time = floor(min_val):1:ceil(max_val)+1;
            bin_centers_time = (bin_edges_time(1:end-1) + bin_edges_time(2:end)) / 2;
            
            [counts_A_time, ~] = histcounts(mRNA_A_time, bin_edges_time, 'Normalization', 'pdf');
            [counts_B_time, ~] = histcounts(mRNA_B_time, bin_edges_time, 'Normalization', 'pdf');
            
            % Plot with increased line width for EPS
            bar(bin_centers_time - bar_width/2, counts_A_time, bar_width, ...
                'FaceColor', color_geneA, 'FaceAlpha', 0.8, ...  % Increased alpha
                'EdgeColor', color_geneA*0.5, 'LineWidth', 1.2);  % Increased line width
            
            bar(bin_centers_time + bar_width/2, counts_B_time, bar_width, ...
                'FaceColor', color_geneB, 'FaceAlpha', 0.8, ...  % Increased alpha
                'EdgeColor', color_geneB*0.5, 'LineWidth', 1.2);  % Increased line width
            
            % Add mean lines with increased width
            y_max_time = max([counts_A_time, counts_B_time]) * 1.1;
            line([mean_A_time, mean_A_time], [0, y_max_time*0.9], ...
                 'Color', color_geneA, 'LineWidth', 2.3, 'LineStyle', '--');  % Increased from 2.0
            line([mean_B_time, mean_B_time], [0, y_max_time*0.9], ...
                 'Color', color_geneB, 'LineWidth', 2.3, 'LineStyle', '--');  % Increased from 2.0
            
            % Formatting
            xlabel('mRNA', 'FontSize', font_size_axis-1, 'FontWeight', 'normal');
            ylabel('P(m)', 'FontSize', font_size_axis-1, 'FontWeight', 'normal');
            title(sprintf('(%s) t = %.1f h (%s)', ...
                  extremal_labels{i}, time_val, extremal_descriptions{i}), ...
                  'FontSize', font_size_title-1, 'FontWeight', 'bold');
            
            % Time point annotation with increased font size
            stat_text_time = {
                sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Gene α:', ...
                        color_geneA(1), color_geneA(2), color_geneA(3))
                sprintf('⟨m⟩ = %.2f', mean_A_time)
                sprintf('σ²/⟨m⟩² = %.4f', noise_A_time)
                ''
                sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Gene β:', ...
                        color_geneB(1), color_geneB(2), color_geneB(3))
                sprintf('⟨m⟩ = %.2f', mean_B_time)
                sprintf('σ²/⟨m⟩² = %.4f', noise_B_time)
            };
            
            text(0.25, 0.95, stat_text_time, ...
                 'Units', 'normalized', 'FontSize', font_size_text, ...  % Original size
                 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
                 'BackgroundColor', 'white', 'Margin', 1, 'Interpreter', 'tex');
            
            grid on;
            set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, ...
                     'Box', 'on', 'LineWidth', 1.2, ...
                     'FontSize', font_size_axis-2);
            xlim([bin_edges_time(1)-0.5, bin_edges_time(end)+0.5]);
        else
            text(0.5, 0.5, 'Data not available', ...
                 'HorizontalAlignment', 'center', 'Units', 'normalized');
        end
    end
    
    %% Optimize figure for EPS output
    % Ensure all text is visible in EPS
    set(findall(fig, '-property', 'FontSize'), 'FontSizeMode', 'manual');
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica');
    
    % Ensure all lines are visible in EPS
    set(findall(fig, 'Type', 'Line'), 'LineWidthMode', 'manual');
    
    %% Save high-quality figures with EPS-specific optimizations
    fprintf('\nSaving high-quality figures...\n');
    
    % Ensure figure is fully rendered
    drawnow;
    
    % Save in multiple formats with EPS-specific settings
    try
        % 1. EPS format for LaTeX - optimized for clarity when enlarged
        % Set up EPS-specific settings
        set(fig, 'Renderer', 'painters');  % Vector renderer for EPS
        set(fig, 'PaperPositionMode', 'auto');
        set(fig, 'InvertHardcopy', 'off');  % Keep exact colors
        
        % Save EPS with specific settings for clarity
        print(fig, '-depsc2', '-r1200', '-painters', 'gillespie_stochastic_simulation.eps');
        fprintf('  EPS saved: gillespie_stochastic_simulation.eps (1200 DPI, vector format)\n');
        fprintf('    Note: EPS is a vector format - will remain clear at any zoom level\n');
        
        % 2. PDF format for publications (vector format)
        print(fig, '-dpdf', '-r1200', '-painters', 'gillespie_stochastic_simulation.pdf');
        fprintf('  PDF saved: gillespie_stochastic_simulation.pdf (1200 DPI, vector format)\n');
        
        % 3. PNG format for general use (high resolution)
        print(fig, '-dpng', '-r600', 'gillespie_stochastic_simulation.png');
        fprintf('  High-res PNG saved: gillespie_stochastic_simulation.png (600 DPI)\n');
        
        % 4. SVG format for editing
        print(fig, '-dsvg', '-r600', 'gillespie_stochastic_simulation.svg');
        fprintf('  SVG saved: gillespie_stochastic_simulation.svg (editable vector)\n');
        
        % 5. FIG format for MATLAB editing
        savefig(fig, 'gillespie_stochastic_simulation.fig');
        fprintf('  FIG saved: gillespie_stochastic_simulation.fig (MATLAB format)\n');
        
        fprintf('\nAll figures saved successfully!\n');
        fprintf('\nFor EPS format (recommended for publications):\n');
        fprintf('  - Uses vector graphics (infinitely scalable)\n');
        fprintf('  - 1200 DPI resolution for high quality\n');
        fprintf('  - All lines and text are vector objects\n');
        fprintf('  - Will remain clear at any magnification\n');
        
    catch ME
        fprintf('Error saving figures: %s\n', ME.message);
        fprintf('Trying alternative save methods...\n');
        
        % Try alternative save methods
        try
            saveas(fig, 'gillespie_stochastic_simulation.eps', 'epsc');
            fprintf('  Saved using saveas: gillespie_stochastic_simulation.eps\n');
        catch
            fprintf('  Could not save figures. Please check write permissions.\n');
        end
    end
    
    %% Display simulation summary
    fprintf('\n=== Simulation Summary ===\n');
    fprintf('Number of simulations: %d\n', num_simulations);
    fprintf('Total simulation time: %.1f hours\n', T_max);
    fprintf('Time step: %.2f hours\n', dt);
    fprintf('\nGene A parameters:\n');
    fprintf('  k_on = %.2f, k_off = %.2f, gamma = %.2f\n', k_on_A, k_off_A, gamma_A);
    fprintf('  V0 = %.2f, a = %.2f, omega = %.2f\n', V0_A, a, omega);
    fprintf('Gene B parameters:\n');
    fprintf('  k_on = %.2f, k_off = %.2f, gamma = %.2f\n', k_on_B, k_off_B, gamma_B);
    fprintf('  k_B = %.2f\n', k_B);
    fprintf('\nGlobal statistics:\n');
    fprintf('  Gene α: mean = %.3f, variance = %.3f, noise = %.5f\n', ...
            mean_A_global, var_A_global, noise_A_global);
    fprintf('  Gene β: mean = %.3f, variance = %.3f, noise = %.5f\n', ...
            mean_B_global, var_B_global, noise_B_global);
    
    %% Tips for EPS usage
    fprintf('\n=== Tips for using EPS in publications ===\n');
    fprintf('1. EPS is a vector format - it uses mathematical formulas to describe shapes\n');
    fprintf('2. Advantages of EPS:\n');
    fprintf('   - Infinitely scalable without loss of quality\n');
    fprintf('   - Small file size for complex figures\n');
    fprintf('   - Perfect for LaTeX documents\n');
    fprintf('   - Professional standard for scientific publications\n');
    fprintf('3. To use in LaTeX: \\includegraphics[width=0.8\\textwidth]{filename.eps}\n');
    fprintf('4. For Word/PowerPoint: Convert EPS to PDF first, then insert\n');
end