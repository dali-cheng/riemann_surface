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
    Re_k = convert_to_0_to_a(angle(z(plot_index)), 2*pi);
    Im_k = -log(abs(z(plot_index)));
    scatter(Re_k / 2 / pi, Im_k, 75, 'k', 'filled'); hold on;
end
hold on;

%% Plot theoretical GBZ
load('./data/gbz_theory.mat'); z = Expression1;
Re_k = convert_to_0_to_a(angle(z), 2*pi);
Im_k = -log(abs(z));
Re_k = [Re_k(511:end); Re_k(1:510)];
Im_k = [Im_k(511:end); Im_k(1:510)];
plot(Re_k / 2 / pi, Im_k, 'LineWidth', 1, 'Color', 'k'); hold on;

% Plot BZ
line([0, 1], [0, 0], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

%% Plot trajectory 2
num_per_roundtrip = 329; % oscilloscope collects 329 datapoints per roundtrip time
pos_1 = round(num_per_roundtrip * 0.6/2);
pos_2 = round(num_per_roundtrip * 1.0/2);

Re_k_traj_1 = transpose(linspace(0.6*pi, pi, pos_2 - pos_1 + 1));
Im_k_traj_1 = 0.2 * ones(1, numel(Re_k_traj_1));
Re_k_traj_2 = transpose(pi * ones(1, 5));
Im_k_traj_2 = [0.3; 0.35; 0.38; 0.39; 0.4];
Re_k_traj_3 = flip(transpose(linspace(0.6*pi, pi, pos_2 - pos_1 + 1)));
Im_k_traj_3 = 0.42 * ones(1, numel(Re_k_traj_3));
Re_k_traj_4 = transpose(0.6*pi * ones(1, 5));
Im_k_traj_4 = [0.4; 0.39; 0.38; 0.35; 0.3];

color_list = hsv(4 * numel(Re_k_traj_1)); color_list = color_list(1 : end/4, :);
scatter(Re_k_traj_1/2/pi, Im_k_traj_1, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_k_traj_1) - 1
    line(Re_k_traj_1(loop_index:loop_index+1)/2/pi, Im_k_traj_1(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_k_traj_1(end), Re_k_traj_2(1)]/2/pi, ...
     [Im_k_traj_1(end), Im_k_traj_2(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

color_list = hsv(4 * numel(Re_k_traj_2)); color_list = color_list(end/4+1 : end/2, :);
scatter(Re_k_traj_2/2/pi, Im_k_traj_2, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_k_traj_2) - 1
    line(Re_k_traj_2(loop_index:loop_index+1)/2/pi, Im_k_traj_2(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_k_traj_2(end), Re_k_traj_3(1)]/2/pi, ...
     [Im_k_traj_2(end), Im_k_traj_3(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

color_list = hsv(4 * numel(Re_k_traj_3)); color_list = color_list(end/2+1 : end/4*3, :);
scatter(Re_k_traj_3/2/pi, Im_k_traj_3, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_k_traj_3) - 1
    line(Re_k_traj_3(loop_index:loop_index+1)/2/pi, Im_k_traj_3(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_k_traj_3(end), Re_k_traj_4(1)]/2/pi, ...
     [Im_k_traj_3(end), Im_k_traj_4(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

color_list = hsv(4 * numel(Re_k_traj_4)); color_list = color_list(end/4*3+1 : end, :);
scatter(Re_k_traj_4/2/pi, Im_k_traj_4, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_k_traj_4) - 1
    line(Re_k_traj_4(loop_index:loop_index+1)/2/pi, Im_k_traj_4(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_k_traj_4(end), Re_k_traj_1(1)]/2/pi, ...
     [Im_k_traj_4(end), Im_k_traj_1(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

%% Plot trajectory 1
pos_1 = round(num_per_roundtrip * 1.0/2);
pos_2 = round(num_per_roundtrip * 1.4/2);

Re_k_traj_1 = linspace(1.0*pi, 1.4*pi, pos_2 - pos_1 + 1);
Im_k_traj_1 = 0.0 * ones(1, numel(Re_k_traj_1));
Re_k_traj_3 = flip(linspace(1.0*pi, 1.4*pi, pos_2 - pos_1 + 1));
Im_k_traj_3 = 0.2 * ones(1, numel(Re_k_traj_3));

Re_k_traj = [Re_k_traj_1, Re_k_traj_3];
Im_k_traj = [Im_k_traj_1, Im_k_traj_3];

color_list = hsv(numel(Re_k_traj));
scatter(Re_k_traj/2/pi, Im_k_traj, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_k_traj) - 1
    line(Re_k_traj(loop_index:loop_index+1)/2/pi, Im_k_traj(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_k_traj(end), Re_k_traj(1)]/2/pi, ...
     [Im_k_traj(end), Im_k_traj(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold off;

%% Figure formats
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
    'YColor', 'k', 'XColor', 'k');
pbaspect([1, 1, 1]);
xlim([0, 1]); ylim([-1, 1] * 0.5);
xticks([0, 1]); yticks([-1, 1] * 0.5);
xticklabels([]); yticklabels([]); 
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');