%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%    Fig4 in the main text   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function main()
    clc; clear all; close all;
    
    % Create figure with publication settings
    fig = figure('Position', [200, 200, 1400, 600], ...
                 'Name', 'Gene Expression Analysis', ...
                 'Color', 'white', ...
                 'Renderer', 'painters', ...  % Vector renderer for high-quality EPS
                 'PaperPositionMode', 'auto', ...
                 'InvertHardcopy', 'off', ...
                 'PaperUnits', 'inches', ...
                 'PaperSize', [14, 6]);  % Larger paper size for better resolution
    
    % Common parameters
    tspan = [0 30];
    y0 = [1; 0; 0; 1; 0; 0; 0; 0; 0; 0];
    
    % Color scheme for different A1 values
    color_A1_0 = [0, 0.5, 0];       % Dark green for A1=0
    color_A1_5 = [0.8, 0, 0.4];     % Magenta for A1=5
    color_A1_10 = [0, 0, 0.5];      % Dark blue for A1=10
    
    % Enhanced line style parameters for EPS clarity
    line_width = 3.0;  % Increased from 2.5 for better visibility in EPS
    solid_style = '-';
    dashed_style = '--';
    
    % ==================== Subplot 1: Average Transcription Levels ====================
    subplot(1, 2, 1);
    hold on; box on; grid on;
    
    % Data storage
    t_data = cell(3, 1);
    mA_data = cell(3, 1);
    mB_data = cell(3, 1);
    
    % Simulation for different A1 values
    [t, Y] = ode45(@myODE_A0_A2_v0, tspan, y0);
    t_data{1} = t; mA_data{1} = Y(:,3); mB_data{1} = Y(:,6);
    
    [t, Y] = ode45(@myODE_A5_A2_v0, tspan, y0);
    t_data{2} = t; mA_data{2} = Y(:,3); mB_data{2} = Y(:,6);
    
    [t, Y] = ode45(@myODE_A10_A2_v0, tspan, y0);
    t_data{3} = t; mA_data{3} = Y(:,3); mB_data{3} = Y(:,6);
    
    % Plotting with increased line widths for EPS
    plot(t_data{1}, mA_data{1}, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', solid_style);
    plot(t_data{2}, mA_data{2}, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', solid_style);
    plot(t_data{3}, mA_data{3}, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', solid_style);
    
    plot(t_data{1}, mB_data{1}, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', dashed_style);
    plot(t_data{2}, mB_data{2}, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', dashed_style);
    plot(t_data{3}, mB_data{3}, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', dashed_style);
    
    % Formatting with increased font sizes for EPS
    xlabel('Time (hours)', 'FontSize', 15, 'FontWeight', 'bold');  % Increased from 14
    ylabel('Mean transcription level', 'FontSize', 15, 'FontWeight', 'bold');  % Increased from 14
    title('(A) Mean Transcription Levels', 'FontSize', 17, 'FontWeight', 'bold');  % Increased from 16
    xlim([0 30]);
    
    % Legend
    legend_handles = gobjects(6, 1);
    legend_handles(1) = plot(NaN, NaN, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', solid_style);
    legend_handles(2) = plot(NaN, NaN, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', solid_style);
    legend_handles(3) = plot(NaN, NaN, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', solid_style);
    legend_handles(4) = plot(NaN, NaN, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', dashed_style);
    legend_handles(5) = plot(NaN, NaN, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', dashed_style);
    legend_handles(6) = plot(NaN, NaN, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', dashed_style);
    
    legend(legend_handles, ...
        {'$m_{\alpha}(t)$, $A_1=0$', ...
         '$m_{\alpha}(t)$, $A_1=5$', ...
         '$m_{\alpha}(t)$, $A_1=10$', ...
         '$m_{\beta}(t)$, $A_1=0$', ...
         '$m_{\beta}(t)$, $A_1=5$', ...
         '$m_{\beta}(t)$, $A_1=10$'}, ...
        'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11, 'NumColumns', 2);  % Increased from 10
    
    ax = gca;
    ax.FontSize = 13;  % Increased from 12
    ax.LineWidth = 1.5;  % Increased from 1.2
    ax.TickDir = 'out';
    grid on; grid minor;
    
    % ==================== Subplot 2: Noise Intensities ====================
    subplot(1, 2, 2);
    hold on; box on;
    
    % Colors for left and right y-axes
    left_ycolor = [0.2, 0.2, 0.2];   % Dark gray for left axis
    right_ycolor = [0.5, 0.2, 0.2];  % Dark red for right axis
    
    % Case 1: A1 = 0
    [t, Y] = ode45(@myODE_A0_A2_v0, tspan, y0);
    mA = Y(:,3); mB = Y(:,6); u1 = Y(:,9); u2 = Y(:,10);
    eta_alpha_1 = (u1 - mA.^2) ./ (mA.^2 + 0.01);
    eta_beta_1 = (u2 - mB.^2) ./ (mB.^2 + 0.01);
    
    yyaxis left;
    h1a = plot(t, eta_alpha_1, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', solid_style);
    ylim([0.5 1]);
    % BOLDED Y-AXIS LABEL - Using \mathbf{} for bold in LaTeX
    ylabel('$\mathbf{CV^2_{\alpha}(t)}$', 'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
    set(gca, 'YColor', left_ycolor, 'FontWeight', 'bold');  % Added FontWeight
    
    yyaxis right;
    h1b = plot(t, eta_beta_1, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', dashed_style);
    ylim([0 3]);
    % BOLDED Y-AXIS LABEL - Using \mathbf{} for bold in LaTeX
    ylabel('$\mathbf{CV^2_{\beta}(t)}$', 'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
    set(gca, 'YColor', right_ycolor, 'FontWeight', 'bold');  % Added FontWeight
    
    % Case 2: A1 = 5
    [t, Y] = ode45(@myODE_A5_A2_v0, tspan, y0);
    mA = Y(:,3); mB = Y(:,6); u1 = Y(:,9); u2 = Y(:,10);
    eta_alpha_2 = (u1 - mA.^2) ./ (mA.^2 + 0.01);
    eta_beta_2 = (u2 - mB.^2) ./ (mB.^2 + 0.01);
    
    yyaxis left;
    h2a = plot(t, eta_alpha_2, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', solid_style);
    
    yyaxis right;
    h2b = plot(t, eta_beta_2, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', dashed_style);
    
    % Case 3: A1 = 10
    [t, Y] = ode45(@myODE_A10_A2_v0, tspan, y0);
    mA = Y(:,3); mB = Y(:,6); u1 = Y(:,9); u2 = Y(:,10);
    eta_alpha_3 = (u1 - mA.^2) ./ (mA.^2 + 0.01);
    eta_beta_3 = (u2 - mB.^2) ./ (mB.^2 + 0.01);
    
    yyaxis left;
    h3a = plot(t, eta_alpha_3, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', solid_style);
    
    yyaxis right;
    h3b = plot(t, eta_beta_3, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', dashed_style);
    
    % Formatting with increased font sizes
    xlabel('Time (hours)', 'FontSize', 15, 'FontWeight', 'bold');  % Increased from 14
    title('(B) Noise Intensities', 'FontSize', 17, 'FontWeight', 'bold');  % Increased from 16
    xlim([0 30]);
    
    % Ensure both y-axis labels are bold (reapply to both sides with explicit bold)
    yyaxis left;
    % Alternative method: Set the ylabel text object properties directly
    ylabel_handle_left = ylabel('$\mathbf{CV^2_{\alpha}(t)}$', 'Interpreter', 'latex', 'FontSize', 15);
    set(ylabel_handle_left, 'FontWeight', 'bold');
    
    yyaxis right;
    ylabel_handle_right = ylabel('$\mathbf{CV^2_{\beta}(t)}$', 'Interpreter', 'latex', 'FontSize', 15);
    set(ylabel_handle_right, 'FontWeight', 'bold');
    
    % Also make axis tick labels bold
    ax = gca;
    ax.FontWeight = 'bold';
    
    set([h1a, h1b, h2a, h2b, h3a, h3b], 'HandleVisibility', 'off');
    
    % ============ MODIFIED: Create legend with ONLY lines (no markers) ============
    % Create dummy lines for legend - WITHOUT markers
    legend_handles2 = gobjects(6, 1);
    legend_handles2(1) = plot(NaN, NaN, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', solid_style, 'Marker', 'none');
    legend_handles2(2) = plot(NaN, NaN, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', solid_style, 'Marker', 'none');
    legend_handles2(3) = plot(NaN, NaN, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', solid_style, 'Marker', 'none');
    legend_handles2(4) = plot(NaN, NaN, 'Color', color_A1_0, 'LineWidth', line_width, 'LineStyle', dashed_style, 'Marker', 'none');
    legend_handles2(5) = plot(NaN, NaN, 'Color', color_A1_5, 'LineWidth', line_width, 'LineStyle', dashed_style, 'Marker', 'none');
    legend_handles2(6) = plot(NaN, NaN, 'Color', color_A1_10, 'LineWidth', line_width, 'LineStyle', dashed_style, 'Marker', 'none');
    
    legend(legend_handles2, ...
        {'$CV^2_{\alpha}(t)$, $A_1=0$', ...
         '$CV^2_{\alpha}(t)$, $A_1=5$', ...
         '$CV^2_{\alpha}(t)$, $A_1=10$', ...
         '$CV^2_{\beta}(t)$, $A_1=0$', ...
         '$CV^2_{\beta}(t)$, $A_1=5$', ...
         '$CV^2_{\beta}(t)$, $A_1=10$'}, ...
        'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11, 'NumColumns', 2);  % Increased from 10
    
    ax = gca;
    ax.FontSize = 13;  % Increased from 12
    ax.LineWidth = 1.5;  % Increased from 1.2
    ax.TickDir = 'out';
    grid on; grid minor;
    
    %% Optimize the figure for EPS output
    % Ensure all text is visible and bold in EPS
    set(findall(fig, '-property', 'FontSize'), 'FontSizeMode', 'manual');
    set(findall(fig, '-property', 'FontWeight'), 'FontWeight', 'bold');
    
    % Ensure all lines are visible in EPS
    set(findall(fig, 'Type', 'Line'), 'LineWidthMode', 'manual');
    
    % Increase axis line widths for better visibility
    all_axes = findall(fig, 'Type', 'axes');
    for i = 1:length(all_axes)
        set(all_axes(i), 'LineWidth', 1.5);  % Increased from 1.2
    end
    
    %% Save high-quality figures with EPS-specific optimizations
    fprintf('Figure generation complete.\n');
    fprintf('\nSaving high-quality figures...\n');
    
    % Ensure figure is fully rendered
    drawnow;
    
    % Save in multiple formats with EPS-specific settings
    try
        % 1. EPS format for LaTeX - optimized for clarity when enlarged
        set(fig, 'Renderer', 'painters');  % Vector renderer for EPS
        set(fig, 'PaperPositionMode', 'auto');
        set(fig, 'InvertHardcopy', 'off');  % Keep exact colors
        
        % Save EPS with high resolution for embedded raster elements
        print(fig, '-depsc2', '-r1200', '-painters', 'gene_expression_analysis.eps');
        fprintf('  EPS saved: gene_expression_analysis.eps (1200 DPI, vector format)\n');
        fprintf('    - Vector format: remains clear at any zoom level\n');
        fprintf('    - All text is bold as specified\n');
        fprintf('    - Line widths increased for better visibility\n');
        fprintf('    - Legend in subplot (B): markers removed, only lines\n');
        
        % 2. PDF format for publications (vector format)
        print(fig, '-dpdf', '-r1200', '-painters', 'gene_expression_analysis.pdf');
        fprintf('  PDF saved: gene_expression_analysis.pdf (1200 DPI, vector format)\n');
        
        % 3. PNG format for general use (high resolution)
        print(fig, '-dpng', '-r600', 'gene_expression_analysis.png');
        fprintf('  High-res PNG saved: gene_expression_analysis.png (600 DPI)\n');
        
        % 4. SVG format for editing
        print(fig, '-dsvg', '-r600', 'gene_expression_analysis.svg');
        fprintf('  SVG saved: gene_expression_analysis.svg (editable vector)\n');
        
        % 5. FIG format for MATLAB editing
        savefig(fig, 'gene_expression_analysis.fig');
        fprintf('  FIG saved: gene_expression_analysis.fig (MATLAB format)\n');
        
        fprintf('\nAll figures saved successfully!\n');
        fprintf('\nEPS Format Features:\n');
        fprintf('  • Vector graphics: Infinite scalability without pixelation\n');
        fprintf('  • All text in subplot (B) y-axis labels is bold\n');
        fprintf('  • Increased line widths (3.0) for better visibility\n');
        fprintf('  • Larger fonts for improved readability\n');
        fprintf('  • Legend in subplot (B): solid and dashed lines only (no markers)\n');
        fprintf('  • Professional quality for academic publications\n');
        
    catch ME
        fprintf('Error saving figures: %s\n', ME.message);
        fprintf('Trying alternative save methods...\n');
        
        % Try alternative save methods
        try
            saveas(fig, 'gene_expression_analysis.eps', 'epsc');
            fprintf('  Saved using saveas: gene_expression_analysis.eps\n');
        catch
            fprintf('  Could not save figures. Please check write permissions.\n');
        end
    end
    
    % Display additional information about the figure
    fprintf('\n=== Figure Details ===\n');
    fprintf('Subplot (A): Mean Transcription Levels\n');
    fprintf('  - Shows m_α(t) and m_β(t) for A1 = 0, 5, 10\n');
    fprintf('  - Solid lines: m_α(t)\n');
    fprintf('  - Dashed lines: m_β(t)\n');
    fprintf('\nSubplot (B): Noise Intensities (CV²)\n');
    fprintf('  - Left y-axis: CV²_α(t) (Noise intensity of gene α)\n');
    fprintf('  - Right y-axis: CV²_β(t) (Noise intensity of gene β)\n');
    fprintf('  - Y-axis labels are bold for emphasis\n');
    fprintf('  - Legend entries use only lines (solid/dashed), no markers\n');
    fprintf('\nParameters: A2 = 2, v2 = 0 for all cases\n');
end

function dYdt = myODE_A0_A2_v0(t, Y)
    % Parameters: A1=0, A2=2, v2=0
    p0 = Y(1); p1 = Y(2); mA = Y(3);
    q0 = Y(4); q1 = Y(5); mB = Y(6);
    m10 = Y(7); m11 = Y(8); u1 = Y(9); u2 = Y(10);
    
    % Gene A parameters
    k_off_A = 0.5;   % Deactivation rate constant for gene α
    k_on_A = 0.5;    % Activation rate constant for gene α
    gamma_A = 1;     % mRNA degradation rate for gene α
    A1 = 0;          % Amplitude of oscillation for gene α
    omega_A = 1;     % Frequency of oscillation
    v1 = 20 + A1 * cos(omega_A * t);  % Transcription rate for gene α
    
    % Gene B parameters
    k_off_B = 0.5;   % Deactivation rate constant for gene β
    k_on_B = 0.5;    % Activation rate constant for gene β 
    gamma_B = 1;     % mRNA degradation rate for gene B
    V0_B = 0;        % Basal transcription rate for gene B
    A2 = 2;          % Coupling strength from gene A to gene B
    
    dp0dt = k_off_A * p1 - k_on_A * p0;
    dp1dt = k_on_A * p0 - k_off_A * p1;
    dmAdt = -gamma_A * mA + v1 * p1;
    dq0dt = k_off_B * q1 - k_on_B * q0;
    dq1dt = k_on_B * p0 - k_off_B * q1;
    dmBdt = -gamma_B * mB + (V0_B + A2 * mA) * q1;
    dm10dt = k_off_A * m11 - (k_on_A + gamma_A) * m10;
    dm11dt = k_on_A * m10 - (k_off_A + gamma_A) * m11 + v1 * p1;
    du1dt = -2 * gamma_A * u1 + v1 * p1 + 2 * v1 * m11 + gamma_A * mA;
    du2dt = (V0_B + A2 * mA) * q1 + ...
            2 * (V0_B^2 + 2 * V0_B * A2 * mA + A2^2 * u1) * ...
            k_off_B * (gamma_B + k_off_B) / ...
            (gamma_B * (k_off_B + k_on_B) * (gamma_B + k_off_B + k_on_B)) + ...
            gamma_B * mB - 2 * gamma_B * u2;
    
    dYdt = [dp0dt; dp1dt; dmAdt; dq0dt; dq1dt; dmBdt; dm10dt; dm11dt; du1dt; du2dt];
end

function dYdt = myODE_A5_A2_v0(t, Y)
    % Parameters: A1=5, A2=2, v2=0
    p0 = Y(1); p1 = Y(2); mA = Y(3);
    q0 = Y(4); q1 = Y(5); mB = Y(6);
    m10 = Y(7); m11 = Y(8); u1 = Y(9); u2 = Y(10);
    
    % Gene A parameters
    k_off_A = 0.5;   
    k_on_A = 0.5;    
    gamma_A = 1;     
    A1 = 5;          
    omega_A = 1;    
    v1 = 20 + A1 * cos(omega_A * t);  
    
    % Gene B parameters
    k_off_B = 0.5;   
    k_on_B = 0.5;    
    gamma_B = 1;     
    V0_B = 0;        
    A2 = 2;          
    
    dp0dt = k_off_A * p1 - k_on_A * p0;
    dp1dt = k_on_A * p0 - k_off_A * p1;
    dmAdt = -gamma_A * mA + v1 * p1;
    dq0dt = k_off_B * q1 - k_on_B * q0;
    dq1dt = k_on_B * p0 - k_off_B * q1;
    dmBdt = -gamma_B * mB + (V0_B + A2 * mA) * q1;
    dm10dt = k_off_A * m11 - (k_on_A + gamma_A) * m10;
    dm11dt = k_on_A * m10 - (k_off_A + gamma_A) * m11 + v1 * p1;
    du1dt = -2 * gamma_A * u1 + v1 * p1 + 2 * v1 * m11 + gamma_A * mA;
    du2dt = (V0_B + A2 * mA) * q1 + ...
            2 * (V0_B^2 + 2 * V0_B * A2 * mA + A2^2 * u1) * ...
            k_off_B * (gamma_B + k_off_B) / ...
            (gamma_B * (k_off_B + k_on_B) * (gamma_B + k_off_B + k_on_B)) + ...
            gamma_B * mB - 2 * gamma_B * u2;
    
    dYdt = [dp0dt; dp1dt; dmAdt; dq0dt; dq1dt; dmBdt; dm10dt; dm11dt; du1dt; du2dt];
end

function dYdt = myODE_A10_A2_v0(t, Y)
    % Parameters: A1=10, A2=2, v2=0
    p0 = Y(1); p1 = Y(2); mA = Y(3);
    q0 = Y(4); q1 = Y(5); mB = Y(6);
    m10 = Y(7); m11 = Y(8); u1 = Y(9); u2 = Y(10);
    
    % Gene A parameters
    k_off_A = 0.5;   
    k_on_A = 0.5;    
    gamma_A = 1;     
    A1 = 10;         
    omega_A = 1;     
    v1 = 20 + A1 * cos(omega_A * t);  
    
    % Gene B parameters
    k_off_B = 0.5;   
    k_on_B = 0.5;    
    gamma_B = 1;    
    V0_B = 0;        
    A2 = 2;          
    
    dp0dt = k_off_A * p1 - k_on_A * p0;
    dp1dt = k_on_A * p0 - k_off_A * p1;
    dmAdt = -gamma_A * mA + v1 * p1;
    dq0dt = k_off_B * q1 - k_on_B * q0;
    dq1dt = k_on_B * p0 - k_off_B * q1;
    dmBdt = -gamma_B * mB + (V0_B + A2 * mA) * q1;
    dm10dt = k_off_A * m11 - (k_on_A + gamma_A) * m10;
    dm11dt = k_on_A * m10 - (k_off_A + gamma_A) * m11 + v1 * p1;
    du1dt = -2 * gamma_A * u1 + v1 * p1 + 2 * v1 * m11 + gamma_A * mA;
    du2dt = (V0_B + A2 * mA) * q1 + ...
            2 * (V0_B^2 + 2 * V0_B * A2 * mA + A2^2 * u1) * ...
            k_off_B * (gamma_B + k_off_B) / ...
            (gamma_B * (k_off_B + k_on_B) * (gamma_B + k_off_B + k_on_B)) + ...
            gamma_B * mB - 2 * gamma_B * u2;
    
    dYdt = [dp0dt; dp1dt; dmAdt; dq0dt; dq1dt; dmBdt; dm10dt; dm11dt; du1dt; du2dt];
end