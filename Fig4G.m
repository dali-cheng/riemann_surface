clearvars; close all; clc;

%% Plot experimental OBC
load('./data/OBC_GBZ.mat', 'imag_k', 'real_E', 'imag_E');
% Load data; values determined manually from the winding intersections

% Apply symmetry
imag_k = [imag_k; -imag_k];
real_E = [real_E; -real_E];
imag_E = [imag_E; -imag_E];

figure; my_colormap = abyss;

for plot_index = 1 : numel(real_E)
    scatter(real_E(plot_index), imag_E(plot_index), 100, ...
        my_colormap(round((imag_k(plot_index) - min(imag_k)) / range(imag_k) * 255 + 1), :), 'filled'); hold on;
end
hold on;

%% Plot theoretical OBC
N = 80; g = 0.2; kappa = 0.1;

H = diag(g * ones(1, N-1), 1) + ...
    diag(g * ones(1, N-1), -1) + ...
    diag(-kappa * ones(1, N-2), 2) + ...
    diag(kappa * ones(1, N-2), -2);
E = eig(H);

% Only take out one branch for plotting
E(imag(E) < 1e-6) = [];
E(real(E) < 0) = [];
E = sort(E, 'descend', 'ComparisonMethod', 'real');

% Manually add in the intersecting point, determined from theoretical winding
intersect_energy = 0.263628;
E(end) = intersect_energy; 

real_E_rescale_factor = 0.04; imag_E_rescale_factor = 0.015;
plot(real(E) * real_E_rescale_factor, imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k'); hold on;
plot(-real(E) * real_E_rescale_factor, imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k'); hold on;
plot(real(E) * real_E_rescale_factor, -imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k'); hold on;
plot(-real(E) * real_E_rescale_factor, -imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k'); hold on;
line([-1, 1] * intersect_energy * real_E_rescale_factor, ...
    [0, 0], [0, 0], 'LineWidth', 1, 'Color', 'k');

%% Figure formats
pbaspect([1, 1, 1]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
    'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
xlim([-1, 1] * 0.02); xticks([-1, 0, 1] * 0.02);
ylim([-1, 1] * 0.004); yticks([-1, 0, 1] * 0.004);
% xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}');
xticklabels([]); yticklabels([]);
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.8], 'Color', 'w');
