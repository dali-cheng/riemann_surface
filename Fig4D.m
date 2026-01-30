clearvars; clc; close all; 
load('./data/sigma_0.39.mat', 'real_E', 'imag_E');
num_per_roundtrip = numel(real_E);
% Extra negative sign here needed because we are using the falling edge of
% the ramp signal in the measurement.
real_E = -real_E; imag_E = -imag_E;

%% Interpolate data using periodic smoothing spline
x = 1 : num_per_roundtrip; p = 0.01; % p is smoothing parameter between 0 and 1
f = spcsp(x, [transpose(real_E); transpose(imag_E)], p);
% This function does periodic smoothing spline.
% https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions
smoothed_data = fnval(f, x);
real_E_smoothed = transpose(smoothed_data(1, :));
imag_E_smoothed = transpose(smoothed_data(2, :));

%% Finding the origin of the real and imaginary parts
imag_E_ctr = mean(imag_E_smoothed);

sigma = 0.39;
b = cosh(sigma) / sinh(2*sigma); % gT_R = 0.2, kappaT_R = 0.1
max_f = 1 + b - 0.5;
if b > 2, min_f = 0.5 - b; else, min_f = (-2 - b^2)/4; end
ratio = -min_f / max_f;
real_E_ctr = min(real_E_smoothed) + ratio/(1+ratio) * range(real_E_smoothed);

real_E = real_E - real_E_ctr;
imag_E = imag_E - imag_E_ctr;
real_E_smoothed = real_E_smoothed - real_E_ctr;
imag_E_smoothed = imag_E_smoothed - imag_E_ctr;

%% Find time origin t = 0
starting_point = convert_to_0_to_a(1 - 0.435 - 65.8 / 360, 1);
% This 65.8 number comes from the GBZ global rotation.
real_E = [real_E; real_E(1 : round(starting_point * num_per_roundtrip))];
real_E(1 : round(starting_point * num_per_roundtrip)) = [];
imag_E = [imag_E; imag_E(1 : round(starting_point * num_per_roundtrip))];
imag_E(1 : round(starting_point * num_per_roundtrip)) = [];
real_E_smoothed = [real_E_smoothed; real_E_smoothed(1 : round(starting_point * num_per_roundtrip))];
real_E_smoothed(1 : round(starting_point * num_per_roundtrip)) = [];
imag_E_smoothed = [imag_E_smoothed; imag_E_smoothed(1 : round(starting_point * num_per_roundtrip))];
imag_E_smoothed(1 : round(starting_point * num_per_roundtrip)) = [];

%% Plot winding
figure;

% Plot the experimental data
scatter(real_E, imag_E, 20, hsv(num_per_roundtrip), 'filled'); hold on;

% Color different sub-regions separately
sub_reg = regions(polyshape(real_E_smoothed, imag_E_smoothed));
plot(sub_reg(1), 'LineWidth', 1, ...
    'EdgeColor', [1, 1, 1] * 0.3, 'FaceColor', [1, 1, 1] * 0.7); hold on;
plot(sub_reg(2), 'LineWidth', 1, ...
    'EdgeColor', [1, 1, 1] * 0.3, 'FaceColor', [1, 1, 1] * 0.7); hold on;
plot(sub_reg(3), 'LineWidth', 1, ...
    'EdgeColor', [1, 1, 1] * 0.3, 'FaceColor', [1, 1, 1] * 0.7); hold on;
plot(sub_reg(4), 'LineWidth', 1, ...
    'EdgeColor', [1, 1, 1] * 0.3, 'FaceColor', [1, 1, 1] * 0.7); hold off;

% Figure formats
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
    'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
pbaspect([1, 1, 1]);
xlim([-1, 1] * 0.004 - 0.0112); ylim([-1, 1] * 0.004);
xticks([-1, 1] * 0.004 - 0.0112); yticks([-1, 1] * 0.004);
xticklabels([]); yticklabels([]); 
% xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');
