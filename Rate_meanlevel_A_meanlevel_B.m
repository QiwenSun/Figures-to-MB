%% Fig2 in main text
%% Three-Curve Plot with Multiple Y-Axes
% This script plots three related curves on the same x-axis with different y-axes
% The curves represent transcription rates and mRNA levels for two genes

clc
clear

%% Time vector
t = 0:0.1:10;  % Time range from 0 to 10 hours with 0.1 hour resolution

%% Gene α parameters (first gene)
V0_A = 20;         % Basal transcription rate for gene α
a = 5;            % Amplitude of periodic transcription
omega = 1;        % Angular frequency of periodic transcription (rad/hour)
k_on_A = 0.5;     % Activation rate constant for gene α
k_off_A = 0.5;    % Deactivation rate constant for gene α
gamma_A = 1;      % mRNA degradation rate for gene α
delay_A = atan(omega/gamma_A);  % Phase delay due to degradation dynamics

%% Gene β parameters (second gene)
V0_B = 0;         % Basal transcription rate for gene β
A2 = 2;           % Proportionality coefficient for gene β synthesis
k_on_B = 0.5;     % Activation rate constant for gene β
k_off_B = 0.5;    % Deactivation rate constant for gene β
gamma_B = 1;      % mRNA degradation rate for gene β
delay_B = atan(omega/gamma_B);  % Phase delay due to degradation dynamics

%% Calculate transcription rate and mRNA levels
% Gene α transcription rate (periodic function)
transcriptionrate_A = V0_A + a*cos(omega*t);

% Gene α mean mRNA level (with delay due to degradation)
meanmRNAlevel_A = (V0_A*k_on_A)/(gamma_A*k_on_A+gamma_A*k_off_A) + ...
    (a*k_on_A)/((k_on_A+k_off_A)*sqrt(gamma_A^2+omega^2))*cos(omega*t - delay_A);

% Gene β mean mRNA level (proportional to gene α with additional delay)
meanmRNAlevel_B = (V0_B + (k_on_B*V0_A*k_on_A)/(gamma_A*k_on_A + gamma_A*k_off_A)) * ...
    (k_on_B/(gamma_B*k_on_B+gamma_B*k_off_B)) + ...
    (a*A2*k_on_A*k_on_B)/((k_on_A+k_off_A)*(k_on_B+k_off_B)* ...
    sqrt((gamma_A^2+omega^2)*(gamma_B^2+omega^2)))*cos(omega*t - delay_A - delay_B);

%% Create figure layout
set(gcf, 'Position', [100, 100, 1200, 500]);  % Set figure size and position
tdl = tiledlayout(1,10);  % Create 1x10 tiled layout for multiple y-axes

% Reduce white space area (can be removed if not needed)
tdl.TileSpacing = 'tight';
tdl.Padding = 'compact';

%% Define color scheme for visualization
color1 = [0, 0.4470, 0.7410];     % Deep blue - gene α transcription rate
color2 = [0.8500, 0.3250, 0.0980]; % Orange-red - gene α mRNA level
color3 = [0.4660, 0.6740, 0.1880]; % Green - gene β mRNA level

%% First axis: Gene α transcription rate ==============================
ax1 = axes(tdl);
hold on
ax1.LineWidth = 1.5;        % Axis line width
ax1.YColor = color1;        % Y-axis color matches curve color
ax1.XLabel.String = 'Time t (hour)';    % X-axis label
ax1.YLabel.String = '$\nu_{\alpha}(t)$'; % Y-axis label: gene α transcription rate
ax1.Layout.TileSpan = [1 9]; % Axis occupies 9/10 width, third y-axis occupies 1/10

% Set label font sizes (increased for better visibility)
ax1.XLabel.FontSize = 16; % X-axis label font size
ax1.YLabel.FontSize = 16; % Y-axis label font size
ax1.YLabel.Color = color1; % Label color matches axis color
ax1.YLabel.Interpreter = 'latex'; % Use LaTeX interpreter for Greek letters

% Set tick font size
ax1.FontSize = 14; % Tick font size

% Set new Y-axis range
ax1.YLim = [13, 27];
ax1.YTick = 13:2:27; % Set tick values

% Plot first data set - gene α transcription rate
plot(ax1, t, transcriptionrate_A, 'Color', color1, 'LineWidth', 2.5);
ylim = get(ax1, 'YLim');

% Add vertical line at x=6.283 (from y-axis minimum to maximum)
line(ax1, [6.283, 6.283], [12, 25], 'Color', color1, 'LineWidth', 1.5, 'LineStyle', '--');
text(ax1, 5.8, ylim(2)-1.5, 'x=6.283', 'Color', color1, 'FontSize', 14, 'FontWeight', 'bold');

%% Second axis: Gene α mRNA level ==============================
ax2 = axes(tdl);
hold on
ax2.LineWidth = 1.5;
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.YColor = color2; % Axis color matches curve color
ax2.YLabel.String = '$m_{\alpha}^{*}(t)$'; % Gene α transcription level
ax2.YLabel.Interpreter = 'latex'; % Use LaTeX interpreter
ax2.Layout.TileSpan = [1 9];

% Plot second data set - gene α mRNA level
plot(ax2, t, meanmRNAlevel_A, 'Color', color2, 'LineWidth', 2.5);
line(ax2, [7.068, 7.068], [7.5, 11.76], 'Color', color2, 'LineWidth', 1.5, 'LineStyle', '--');
text(ax2, 6.7, 11.98, 'x=7.068', 'Color', color2, 'FontSize', 14, 'FontWeight', 'bold');

% Set label font sizes
ax2.YLabel.FontSize = 16; % Y-axis label font size
ax2.YLabel.Color = color2; % Label color matches axis color
ax2.FontSize = 14; % Tick font size

% Set new Y-axis range
ax2.YLim = [7.5, 12.5];
ax2.YTick = 7.5:0.5:12.5; % Set tick values

%% Third axis: Gene β mRNA level ==============================
ax3 = axes(tdl);
hold on
ax3.Color = 'none';
ax3.YColor = 'none';
ax3.Layout.TileSpan = [1 9];

% Plot third data set - gene β mRNA level
plot(ax3, t, meanmRNAlevel_B, 'Color', color3, 'LineWidth', 2.5);
line(ax3, [7.853, 7.853], [0.5, 3.75], 'Color', color3, 'LineWidth', 1.5, 'LineStyle', '--');
text(ax3, 7.6, 3.92, 'x=7.853', 'Color', color3, 'FontSize', 14, 'FontWeight', 'bold');

% Set tick font size
ax3.FontSize = 14; % Tick font size

% Set new Y-axis range
ax3.YLim = [3.5, 10];
ax3.YTick = 4:2:10; % Set tick values

% Link x-axes of all coordinate areas
linkaxes(tdl.Children,'x')

%% Create third y-axis for gene β mRNA level (occupying 1/10 width)
ax4 = axes(tdl, 'LineWidth', 1.5, 'YAxisLocation', 'right', ...
    'Color', 'none', 'XColor', 'none');
ax4.YColor = color3; % Axis color matches curve color
ax4.YLabel.String = '$m_{\beta}^{*}(t)$'; % Gene β transcription level
ax4.YLabel.Interpreter = 'latex'; % Use LaTeX interpreter
ax4.Layout.Tile = 'east';
ax4.Layout.Tile = 10;
linkaxes([ax3,ax4],'y')

% Set label font sizes
ax4.YLabel.FontSize = 16; % Y-axis label font size
ax4.YLabel.Color = color3; % Label color matches axis color
ax4.FontSize = 14; % Tick font size

% Set new Y-axis range
ax4.YLim = [0.7, 4.3];
ax4.YTick = 1:0.5:4; % Set tick values

%% Set x-axis limits and labels
xlim([0, 10]);
xlabel(ax1, 'Time t (hours)', 'FontSize', 16);

%% Create legend
% Create dummy graphic objects for legend
h1 = plot(ax1, NaN, NaN, 'Color', color1, 'LineWidth', 2.5);
h2 = plot(ax2, NaN, NaN, 'Color', color2, 'LineWidth', 2.5);
h3 = plot(ax3, NaN, NaN, 'Color', color3, 'LineWidth', 2.5);

% Create legend with LaTeX syntax
lgd = legend([h1, h2, h3], {'$\nu_{\alpha}(t)$', '$m_{\alpha}^{*}(t)$', '$m_{\beta}^{*}(t)$'}, ...
    'Location', 'northeast', 'FontSize', 14, 'Box', 'off', 'Interpreter', 'latex');

%% Add grid to all axes
grid(ax1, 'on');
grid(ax2, 'on');
grid(ax3, 'on');

% Set grid style
set([ax1, ax2, ax3], 'GridLineStyle', ':', 'GridAlpha', 0.3);

%% Beautify figure
set(gcf, 'Color', 'w'); % White background

% Add title
%title(ax1, 'Transcription Rates and mRNA Levels for Genes α and β', ...
 %   'FontSize', 16, 'FontWeight', 'bold');

%% Save high-quality figure
fprintf('Saving high-quality figures...\n');

try
    % 1. PNG format for general use (high resolution)
    print('-dpng', '-r600', 'fig2_transcription_rates_mRNA_levels.png');
    fprintf('  High-res PNG saved: fig2_transcription_rates_mRNA_levels.png (600 DPI)\n');
    
    % 2. PDF format for publications (vector format)
    print('-dpdf', '-r600', 'fig2_transcription_rates_mRNA_levels.pdf');
    fprintf('  PDF saved: fig2_transcription_rates_mRNA_levels.pdf (vector format)\n');
    
    % 3. SVG format for editing
    print('-dsvg', '-r600', 'fig2_transcription_rates_mRNA_levels.svg');
    fprintf('  SVG saved: fig2_transcription_rates_mRNA_levels.svg (editable vector)\n');
    
    % 4. EPS format for LaTeX
    print('-depsc2', '-r600', 'fig2_transcription_rates_mRNA_levels.eps');
    fprintf('  EPS saved: fig2_transcription_rates_mRNA_levels.eps (for LaTeX)\n');
    
    fprintf('\nFigure 2 saved successfully!\n');
    
catch ME
    fprintf('Error saving figures: %s\n', ME.message);
    fprintf('Trying alternative save methods...\n');
    
    % Try alternative save methods
    try
        saveas(gcf, 'fig2_transcription_rates_mRNA_levels.png');
        fprintf('  Saved using saveas: fig2_transcription_rates_mRNA_levels.png\n');
    catch
        fprintf('  Could not save figures. Please check write permissions.\n');
    end
end

%% Display parameter summary
fprintf('\n=== Parameter Summary for Figure 2 ===\n');
fprintf('Time range: 0 to 10 hours, resolution: 0.1 hour\n');
fprintf('\nGene α parameters:\n');
fprintf('  Basal transcription rate (V0_A): %.1f\n', V0_A);
fprintf('  Amplitude (a): %.1f\n', a);
fprintf('  Angular frequency (ω): %.1f rad/hour\n', omega);
fprintf('  Activation rate (k_on_A): %.1f\n', k_on_A);
fprintf('  Deactivation rate (k_off_A): %.1f\n', k_off_A);
fprintf('  mRNA degradation rate (γ_A): %.1f\n', gamma_A);
fprintf('  Phase delay (δ_A): %.4f radians\n', delay_A);
fprintf('\nGene β parameters:\n');
fprintf('  Basal transcription rate (V0_B): %.1f\n', V0_B);
fprintf('  Proportionality coefficient (A2): %.1f\n', A2);
fprintf('  Activation rate (k_on_B): %.1f\n', k_on_B);
fprintf('  Deactivation rate (k_off_B): %.1f\n', k_off_B);
fprintf('  mRNA degradation rate (γ_B): %.1f\n', gamma_B);
fprintf('  Phase delay (δ_B): %.4f radians\n', delay_B);
fprintf('  Total phase delay (δ_A + δ_B): %.4f radians\n', delay_A + delay_B);
fprintf('\nVertical line positions:\n');
fprintf('  ν_α(t) peak: x = 6.283 (2π)\n');
fprintf('  m_α*(t) peak: x = 7.068 (2π + δ_A)\n');
fprintf('  m_β*(t) peak: x = 7.853 (2π + δ_A + δ_B)\n');