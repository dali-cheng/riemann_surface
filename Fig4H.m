clearvars; close all; clc;

%% Plot experimental GBZ
load('./data/OBC_GBZ.mat', 'imag_k', 'real_k_1', 'real_k_2');
% Load data; values determined manually from the winding intersections

rot_angle = 65.8; % Global time delay from optical signal vs electrical signal

% Apply symmetry
real_k = [real_k_1; real_k_2];
imag_k = [imag_k; imag_k];

z = exp(1i * (real_k + 1i * imag_k)) * exp(1i * deg2rad(rot_angle));
z = [z; -1./z];
imag_k = [imag_k; -imag_k];

figure; my_colormap = abyss;

for plot_index = 1 : numel(z)
    scatter(real(z(plot_index)), imag(z(plot_index)), 100, ...
        my_colormap(round((imag_k(plot_index) - min(imag_k)) / range(imag_k) * 255 + 1), :), 'filled'); hold on;
end
hold on;

%% Plot theoretical GBZ
load('./data/gbz_theory.mat'); z = Expression1;
plot(real(z), imag(z), 'LineWidth', 1, 'Color', 'k'); hold on;
line([real(z(1)), real(z(end))], [imag(z(1)), imag(z(end))], ...
    'LineWidth', 1, 'Color', 'k');

%% Figure formats
% Plot unit circle
rectangle('Position', [-1, -1, 2, 2], 'Curvature', [1, 1], ...
    'LineWidth', 1, 'LineStyle', '--');

axis equal;
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off');
xlim([-1, 1] * 2); xticks([-2, 0, 2]);
ylim([-1, 1] * 2); yticks([-2, 0, 2]);
% xlabel('Re({\itz})'); ylabel('Im({\itz})');
xticklabels([]); yticklabels([]);
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.8], 'Color', 'w');


exportgraphics(gcf, 'Figure4H.pdf', 'ContentType', 'vector');