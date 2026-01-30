clearvars; clc; close all; 
load('./data/sigma_0.39.mat');

%% Find time origin t = 0
starting_point = convert_to_0_to_a(1 - 0.435 - 65.8 / 360, 1);
% This 65.8 number comes from the GBZ global rotation
num_per_roundtrip = size(bandstr, 2);
bandstr = [bandstr, bandstr(:, 1 : round(starting_point * num_per_roundtrip))];
bandstr(:, 1 : round(starting_point * num_per_roundtrip)) = [];

%% Interpolate data using periodic smoothing spline
x = 1 : numel(real_E); p = 0.01; % p is smoothing parameter between 0 and 1
f = spcsp(x, [transpose(real_E); transpose(imag_E)], p);
% This function does periodic smoothing spline.
% https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions
smoothed_data = fnval(f, x);
real_E_smoothed = transpose(smoothed_data(1, :));
imag_E_smoothed = transpose(smoothed_data(2, :));

%% Finding the origin of the real part
sigma = 0.39;
b = cosh(sigma) / sinh(2*sigma); % gT_R = 0.2, kappaT_R = 0.1
max_f = 1 + b - 0.5;
if b > 2, min_f = 0.5 - b; else, min_f = (-2 - b^2)/4; end
ratio = -min_f / max_f;
real_E_ctr = min(real_E_smoothed) + ratio/(1+ratio) * range(real_E_smoothed);

%% Generaing plotting axes
fast_time = linspace(0, 1, num_per_roundtrip);
delta_omega = (0 : (size(bandstr, 1) - 1)) / delta_x_mean - real_E_ctr;
delta_omega = -delta_omega;
% We need to add a negative sign to the delta_omega, because
% we record data on the falling edge of the ramp signal.
[fast_time_plot, delta_omega_plot] = meshgrid(fast_time, delta_omega);

%% Plot bandstr
figure;
surf(fast_time_plot, delta_omega_plot, bandstr, 'EdgeColor', 'none');

xlim([0, 1]); ylim([-1, 1] * 0.2); clim([0.162, 0.181]); grid off;
xticks([0, 0.25, 0.5, 0.75, 1]); yticks([-1, -0.5, 0, 0.5, 1] * 0.2);
% xlabel('{\itt} / {\itT_R}'); ylabel('\delta{\it\omega} / {\it\Omega_R}');
xticklabels([]); yticklabels([]);
colormap(flipud(jet)); view([0, 0, 1]); 
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'TickDir', 'out', 'labelfontsizemultiplier', 1, 'linewidth', 1, ...
    'Layer', 'Top', 'Box', 'off');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');