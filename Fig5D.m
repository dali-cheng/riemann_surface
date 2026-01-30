clearvars; clc; close all;
num_per_roundtrip = 329; % oscilloscope collects 329 datapoints per roundtrip time
pos_1 = round(num_per_roundtrip * 0.6/2);
pos_2 = round(num_per_roundtrip * 1.0/2);

%% Rectangular section 1: Im(k) = 0.2, Re(k) from 0.6pi to pi
[real_E, imag_E] = data_process(0.2, './data/sigma_0.2.mat', 0.2);
Re_E_traj_1 = real_E(pos_1 : pos_2);
Im_E_traj_1 = imag_E(pos_1 : pos_2);

%% Rectangular section 3: Im(k) = 0.42, Re(k) from pi to 0.6pi
[real_E, imag_E] = data_process(0.42, './data/sigma_0.42.mat', 0.615);
Re_E_traj_3 = flip(real_E(pos_1 : pos_2));
Im_E_traj_3 = flip(imag_E(pos_1 : pos_2));

%% Rectangular section 2 or 4: Re(k) = 0.6pi or pi, Im(k) in {0.3, 0.35, 0.38, 0.39, 0.4}
Re_E_traj_2 = []; Im_E_traj_2 = []; Re_E_traj_4 = []; Im_E_traj_4 = [];

% sigma = 0.3
[real_E, imag_E] = data_process(0.3, './data/sigma_0.3.mat', 0.475);
Re_E_traj_2 = [Re_E_traj_2; real_E(pos_2)]; Im_E_traj_2 = [Im_E_traj_2; imag_E(pos_2)];
Re_E_traj_4 = [Re_E_traj_4; real_E(pos_1)]; Im_E_traj_4 = [Im_E_traj_4; imag_E(pos_1)];

% sigma = 0.35
[real_E, imag_E] = data_process(0.35, './data/sigma_0.35.mat', 0.3);
Re_E_traj_2 = [Re_E_traj_2; real_E(pos_2)]; Im_E_traj_2 = [Im_E_traj_2; imag_E(pos_2)];
Re_E_traj_4 = [Re_E_traj_4; real_E(pos_1)]; Im_E_traj_4 = [Im_E_traj_4; imag_E(pos_1)];

% sigma = 0.38
[real_E, imag_E] = data_process(0.38, './data/sigma_0.38.mat', 0.805);
Re_E_traj_2 = [Re_E_traj_2; real_E(pos_2)]; Im_E_traj_2 = [Im_E_traj_2; imag_E(pos_2)];
Re_E_traj_4 = [Re_E_traj_4; real_E(pos_1)]; Im_E_traj_4 = [Im_E_traj_4; imag_E(pos_1)];

% sigma = 0.39
[real_E, imag_E] = data_process(0.39, './data/sigma_0.39.mat', 0.435);
Re_E_traj_2 = [Re_E_traj_2; real_E(pos_2)]; Im_E_traj_2 = [Im_E_traj_2; imag_E(pos_2)];
Re_E_traj_4 = [Re_E_traj_4; real_E(pos_1)]; Im_E_traj_4 = [Im_E_traj_4; imag_E(pos_1)];

% sigma = 0.4
[real_E, imag_E] = data_process(0.4, './data/sigma_0.4.mat', 0.12);
Re_E_traj_2 = [Re_E_traj_2; real_E(pos_2)]; Im_E_traj_2 = [Im_E_traj_2; imag_E(pos_2)];
Re_E_traj_4 = [Re_E_traj_4; real_E(pos_1)]; Im_E_traj_4 = [Im_E_traj_4; imag_E(pos_1)];

Re_E_traj_4 = flip(Re_E_traj_4); Im_E_traj_4 = flip(Im_E_traj_4);

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
    'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); hold on;
plot(-real(E) * real_E_rescale_factor, imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); hold on;
plot(real(E) * real_E_rescale_factor, -imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); hold on;
plot(-real(E) * real_E_rescale_factor, -imag(E) * imag_E_rescale_factor, ...
    'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); hold on;
line([-1, 1] * intersect_energy * real_E_rescale_factor, ...
    [0, 0], [0, 0], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); hold on;

%% Plot encircling EP
color_list = hsv(4 * numel(Re_E_traj_1)); color_list = color_list(1 : end/4, :);
scatter(Re_E_traj_1, Im_E_traj_1, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_E_traj_1) - 1
    line(Re_E_traj_1(loop_index:loop_index+1), Im_E_traj_1(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_E_traj_1(end), Re_E_traj_2(1)], ...
     [Im_E_traj_1(end), Im_E_traj_2(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

color_list = hsv(4 * numel(Re_E_traj_2)); color_list = color_list(end/4+1 : end/2, :);
scatter(Re_E_traj_2, Im_E_traj_2, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_E_traj_2) - 1
    line(Re_E_traj_2(loop_index:loop_index+1), Im_E_traj_2(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_E_traj_2(end), Re_E_traj_3(1)], ...
     [Im_E_traj_2(end), Im_E_traj_3(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

color_list = hsv(4 * numel(Re_E_traj_3)); color_list = color_list(end/2+1 : end/4*3, :);
scatter(Re_E_traj_3, Im_E_traj_3, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_E_traj_3) - 1
    line(Re_E_traj_3(loop_index:loop_index+1), Im_E_traj_3(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_E_traj_3(end), Re_E_traj_4(1)], ...
     [Im_E_traj_3(end), Im_E_traj_4(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

color_list = hsv(4 * numel(Re_E_traj_4)); color_list = color_list(end/4*3+1 : end, :);
scatter(Re_E_traj_4, Im_E_traj_4, 20, color_list, "filled"); hold on;
for loop_index = 1 : numel(Re_E_traj_4) - 1
    line(Re_E_traj_4(loop_index:loop_index+1), Im_E_traj_4(loop_index:loop_index+1), ...
        'LineWidth', 1, 'Color', color_list(loop_index, :)); hold on;
end
line([Re_E_traj_4(end), Re_E_traj_1(1)], ...
     [Im_E_traj_4(end), Im_E_traj_1(1)], ...
        'LineWidth', 1, 'Color', color_list(end, :)); hold on;

%% Figure formats
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
    'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
pbaspect([1, 1, 1]);
xlim([-1, 1] * 0.004 - 0.0112); ylim([-1, 1] * 0.004);
xticks([-1, 1] * 0.004 - 0.0112); yticks([-1, 1] * 0.004);
xticklabels([]); yticklabels([]); 
% xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');


function [real_E, imag_E] = data_process(sigma, filename, stpt)
load(filename, 'real_E', 'imag_E');
num_per_roundtrip = numel(real_E); real_E = -real_E; imag_E = -imag_E;
% Extra negative sign here needed because we are using the falling edge of
% the ramp signal in the measurement.

% Interpolate data using periodic smoothing spline
x = 1 : num_per_roundtrip; p = 0.01; % p is smoothing parameter between 0 and 1
f = spcsp(x, [transpose(real_E); transpose(imag_E)], p);
% This function does periodic smoothing spline.
% https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions
smoothed_data = fnval(f, x);
real_E_smoothed = transpose(smoothed_data(1, :));
imag_E_smoothed = transpose(smoothed_data(2, :));

% Finding the origin of the real and imaginary parts
imag_E_ctr = mean(imag_E_smoothed);
imag_E = imag_E - imag_E_ctr;

if sigma == 0, ratio = 1;
else
    b = cosh(sigma) / sinh(2 * sigma);
    max_f = 1 + b - 0.5;
    if b > 2, min_f = 0.5 - b; else, min_f = (-2 - b^2)/4; end
    ratio = -min_f / max_f;
end
real_E_ctr = min(real_E_smoothed) + ratio/(1+ratio) * range(real_E_smoothed);
real_E = real_E - real_E_ctr;

% Find time origin t = 0
starting_point = convert_to_0_to_a(1 - stpt - 65.8 / 360, 1);
% This 65.8 number comes from the GBZ global rotation.
real_E = [real_E; real_E(1 : round(starting_point * num_per_roundtrip))];
real_E(1 : round(starting_point * num_per_roundtrip)) = [];
imag_E = [imag_E; imag_E(1 : round(starting_point * num_per_roundtrip))];
imag_E(1 : round(starting_point * num_per_roundtrip)) = [];
end