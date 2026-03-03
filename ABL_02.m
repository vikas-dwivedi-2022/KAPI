clc; clear; close all

nu = 0.001;

% List of Nc values
Nc_list = [500, 1000, 2000, 3000, 4000, 5000];

% Marker styles
marker_styles = {'-o', '-s', '-^', '-d', '-v', '-p'};

% Create wide figure
f = figure('Units','inches','Position',[1 1 9 5]);
hold on; grid on; box on;

% Preallocate line handles
h = gobjects(length(Nc_list),1);

% Plot curves
for i = 1:length(Nc_list)
    Nc = Nc_list(i);
    filename = sprintf('logFactor_logJ_Nc%d.csv', Nc);

    if ~isfile(filename)
        warning('File %s not found. Skipping.', filename);
        continue
    end

    data = readmatrix(filename);
    logFactor = data(:,1);
    logJ = data(:,2);

    h(i) = plot(logFactor, logJ, marker_styles{i}, ...
        'LineWidth', 2, 'MarkerSize', 6);
end

% Axis labels (use single $ for MATLAB LaTeX)
xlabel('$\log_{10}(\sigma_F)$','Interpreter','latex','FontSize',24);
ylabel('$\log_{10}(J)$','Interpreter','latex','FontSize',24);

% Title (clean LaTeX-safe formatting)
title(sprintf('Baseline PI-ELM, $\\nu = %.3g$', nu), ...
      'Interpreter','latex','FontSize',20);
% Create legend strings
legend_strings = arrayfun(@(Nc) sprintf('$N_{rbf} = %d$', Nc), ...
                          Nc_list, 'UniformOutput', false);

% Add legend using handles (robust method)
lgd = legend(h, legend_strings, ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'Location','eastoutside');

% Improve appearance
set(gca, 'FontSize', 18, ...
         'TickLabelInterpreter','latex', ...
         'LineWidth',1.2);

% Export high-resolution figure
exportgraphics(f, 'FIG_J_PIELM_Vs_eta_Vs_N-RBF.png', 'Resolution', 300);