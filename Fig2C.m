clearvars; close all; clc;
g = 0.2; kappa = 0.1;

%% Calculate Riemann surface
real_k = linspace(0, 2, 201) * pi;
imag_k = linspace(-1, 1, 201) * 0.7;
[real_k_plot, imag_k_plot] = meshgrid(real_k, imag_k);
k = real_k_plot + 1i * imag_k_plot;

E = g * exp(-1i*k) + g * exp(1i*k) + kappa * exp(-2i*k) - kappa * exp(2i*k);

figure; 
my_colormap = hsv(numel(real_k));
my_color = zeros(size(E, 1), size(E, 2), 3);
for idx1 = 1 : size(E, 1)
    for idx2 = 1 : size(E, 2)
        [~, loc] = find_closest(real_k, real_k_plot(idx1, idx2));
        my_color(idx1, idx2, :) = my_colormap(loc, :);
    end
end
surf(real(E), imag(E), imag_k_plot, ...
    my_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;

%% Figure formats
lightangle(180, 5); lighting none; grid off;
% xlabel('Re({\itE})'); ylabel('Im({\itE})'); zlabel('Im({\itk})'); 
xticklabels([]); yticklabels([]); zticklabels([]);
xlim([-1, 1] * 1); ylim([-1, 1] * 1); zlim([-1, 1] * 0.7); view([-38.2, 11]);
set(gca, 'fontname', 'Arial', 'fontsize', 28, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'box', 'off', ...
    'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
xticks([-1, 0, 1] * 1); yticks([-1, 0, 1] * 1); zticks([-1, -0.5, 0, 0.5, 1] * 0.7);
set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.45 0.8], 'color', 'w');