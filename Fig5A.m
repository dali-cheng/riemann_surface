clearvars; clc; close all; 
load('./data/sigma_0.38.mat', 'real_E', 'imag_E');
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

sigma = 0.38;
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
starting_point = convert_to_0_to_a(1 - 0.805 - 65.8 / 360, 1);
% This 65.8 number comes from the GBZ global rotation.
real_E = [real_E; real_E(1 : round(starting_point * num_per_roundtrip))];
real_E(1 : round(starting_point * num_per_roundtrip)) = [];
imag_E = [imag_E; imag_E(1 : round(starting_point * num_per_roundtrip))];
imag_E(1 : round(starting_point * num_per_roundtrip)) = [];
real_E_smoothed = [real_E_smoothed; real_E_smoothed(1 : round(starting_point * num_per_roundtrip))];
real_E_smoothed(1 : round(starting_point * num_per_roundtrip)) = [];
imag_E_smoothed = [imag_E_smoothed; imag_E_smoothed(1 : round(starting_point * num_per_roundtrip))];
imag_E_smoothed(1 : round(starting_point * num_per_roundtrip)) = [];

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

figure;
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

%% Plot winding
scatter(real_E, imag_E, 20, hsv(num_per_roundtrip), 'filled'); hold off;

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
