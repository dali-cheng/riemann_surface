clearvars; clc; close all;
sigma_list = [0.0, 0.1, 0.2, 0.3, 0.35, 0.38, 0.39, ...
    0.4, 0.42, 0.45, 0.5, 0.55, 0.6, 0.7];
time_origin_list = convert_to_0_to_a([0.69, 0.03, 0.2, 0.475, 0.3, ...
    0.805, 0.435, 0.12, 0.615, 0.57, 0.84, 0.792, 0.492, 0.09] + 65.8 / 360, 1);
% This 65.8 number comes from the GBZ global rotation.

figure;
for sigma_index = 1 : numel(sigma_list)
    sigma = sigma_list(sigma_index);

    %% Load data
    load(['./data/sigma_', num2str(sigma), '.mat'], 'real_E', 'imag_E');
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

    if sigma == 0, ratio = 1;
    else
        b = cosh(sigma) / sinh(2 * sigma);
        max_f = 1 + b - 0.5;
        if b > 2, min_f = 0.5 - b; else, min_f = (-2 - b^2)/4; end
        ratio = -min_f / max_f;
    end
    real_E_ctr = min(real_E_smoothed) + ratio/(1+ratio) * range(real_E_smoothed);

    real_E = real_E - real_E_ctr;
    imag_E = imag_E - imag_E_ctr;
    real_E_smoothed = real_E_smoothed - real_E_ctr;
    imag_E_smoothed = imag_E_smoothed - imag_E_ctr;

    %% Find time origin t = 0
    starting_point = 1 - time_origin_list(sigma_index);

    real_E = [real_E; real_E(1 : round(starting_point * num_per_roundtrip))];
    real_E(1 : round(starting_point * num_per_roundtrip)) = [];
    imag_E = [imag_E; imag_E(1 : round(starting_point * num_per_roundtrip))];
    imag_E(1 : round(starting_point * num_per_roundtrip)) = [];
    real_E_smoothed = [real_E_smoothed; real_E_smoothed(1 : round(starting_point * num_per_roundtrip))];
    real_E_smoothed(1 : round(starting_point * num_per_roundtrip)) = [];
    imag_E_smoothed = [imag_E_smoothed; imag_E_smoothed(1 : round(starting_point * num_per_roundtrip))];
    imag_E_smoothed(1 : round(starting_point * num_per_roundtrip)) = [];

    %% Plot winding
    color_list = hsv(num_per_roundtrip);
    scatter3(real_E, imag_E, sigma * ones(size(real_E)), 15, color_list, 'filled'); hold on;

    %% Plot symmetrical winding
    % E -> -E, sigma -> -sigma, Re(k) -> pi - Re(k)
    if sigma > 0
        color_list = [color_list; color_list(1 : round(num_per_roundtrip/2), :)];
        color_list(1 : round(num_per_roundtrip/2), :) = [];
        scatter3(-real_E, -imag_E, -sigma * ones(size(real_E)), 15, ...
            flipud(color_list), 'filled'); hold on;
    end
end

constantplane([0, 0, 1], 0, ...
    'FaceColor', [1, 1, 1] * 0.75, 'FaceAlpha', 0.6); hold on;
constantplane([0, 0, 1], 0.39, ...
    'FaceColor', [1, 1, 1] * 0.75, 'FaceAlpha', 0.6); hold on;

%% Figure formats
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'box', 'off', ...
    'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
xlim([-1, 1] * 0.04);  xticks([-1, 0, 1] * 0.04);
ylim([-1, 1] * 0.015); yticks([-1, 0, 1] * 0.015);
zlim([-1, 1] * 0.7);   zticks([-1, -0.5, 0, 0.5, 1] * 0.7);
xticklabels([]); yticklabels([]); zticklabels([]);
% xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}'); zlabel('Im({\itk})');
grid off; view([-38.2, 11]);
set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.45 0.8], 'color', 'w');