clearvars; clc; close all; 
load('./data/sigma_0.mat');

%% Find time origin t = 0
starting_point = convert_to_0_to_a(1 - 0.69 - 65.8 / 360, 1);
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
real_E_ctr = (min(real_E_smoothed) + max(real_E_smoothed)) / 2;
% This only applies to sigma = 0

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