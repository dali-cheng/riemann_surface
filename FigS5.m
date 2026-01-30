clearvars; clc; close all; 

%% Definition of parameters used in the simulation

gamma_coupling = 0.05; 
% Power splitting ratio of fiber coupler

g_TR = 0.045; kappa_TR = 0.02; 
% We chose these values for optimal agreement with experimental results.
% The ratio g/kappa is not exactly 2:1, and the slight deviation could be
% ascribed to imperfect electronics and/or modulators, where the AWG 
% voltages could translate to different modulation strengths.

gamma0 = 0.5; 
% Round-trip field loss of the resonator

sigma_list = [
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, ...
    0.35, 0.36, 0.37, 0.38, 0.39, 0.4, ...
    0.45, 0.5, 0.55, 0.6, 0.65, 0.7
    ]; % Set of Im(k) values to be probed

delta_omega_min = -0.2;
delta_omega_max = 1.2;
% Range of interest of the laser detuning (unit: FSR)

num_pts_within_roundtrip = 350; num_roundtrips = 6401;
% How many datapoints does each round-trip time contain, and
% how many round-trips do we simulate

bandstr = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(bandstr) + 1) * 2*pi; 
Omegat_list(end) = [];
delta_omega_list = linspace(delta_omega_min, delta_omega_max, numel(bandstr) + 1) * 2*pi; 
delta_omega_list(end) = [];

riemann_surface_real_E = dictionary;
riemann_surface_imag_E = dictionary;
% Storage for Re(E) and Im(E)

for sigma = sigma_list

%% Simulate transmission spectrum
    for loop_index = 1 : numel(bandstr)
        % Calculate transmission at each time point
        Omegat = Omegat_list(loop_index);
        delta_omega = delta_omega_list(loop_index);

        PM = 2 * g_TR * cosh(sigma) * cos(Omegat) ...
            + 2 * kappa_TR * sinh(2*sigma) * cos(2*Omegat);
        AM = 2 * g_TR * sinh(sigma) * sin(Omegat) ...
            + 2 * kappa_TR * cosh(2*sigma) * sin(2*Omegat);
        T = exp(-1i * PM) * exp(-AM);
        P = exp(1i * delta_omega) * exp(-gamma0);
    
        psi_out = sqrt(1 - gamma_coupling) ...
            - gamma_coupling * T * P / (1 - T * P * sqrt(1 - gamma_coupling));
        % Steady-state resonator transmittance based on input-output
        % coupling formalism
    
        bandstr(loop_index) = abs(psi_out)^2;
    end
    bandstr = transpose(reshape(bandstr, num_pts_within_roundtrip, num_roundtrips));
    bandstr = bandstr / max(max(bandstr)); % Normalize

    %% Plot Fig. S5(A) and S5(C)
    if (sigma == 0.0) || (sigma == 0.39)
        figure;
        
        [Omegat_plot, delta_omega_plot] = meshgrid( ...
            linspace(0, 1, num_pts_within_roundtrip), ...
            linspace(delta_omega_min, delta_omega_max, num_roundtrips) ...
            );
        surf(Omegat_plot, delta_omega_plot, bandstr, 'EdgeColor', 'none');
        
        xlim([0, 1]); ylim([-1, 1] * 0.2); grid off;
        clim([0.78, 1]);
        xticks([0, 0.25, 0.5, 0.75, 1]); yticks([-1, -0.5, 0, 0.5, 1] * 0.2);
        % xlabel('{\itt} / {\itT_R}'); ylabel('\delta{\it\omega} / {\it\Omega_R}');
        xticklabels([]); yticklabels([]);
        colormap(flipud(jet)); view([0, 0, 1]); 
        set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
            'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off');
        set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');
    end

    %% Fitting
    % Experimental results are obtained using the same Lorentzian fitting
    num_res = 2; % The number of resonances inside the delta_omega range of interest
    res_loc = zeros(num_res, num_pts_within_roundtrip);
    res_wid = zeros(1, num_pts_within_roundtrip);
    delta_x_mean = 0; % Will provide a more accurate conversion from matrix index to FSR
    
    for col_index = 1 : num_pts_within_roundtrip
        % Send column slices to fitting; 
        % xi(delta_omega) should have Lorentzian lineshape
        if col_index == 1
            % Fit with default initial value
            [res_loc(:, col_index), res_wid(col_index), ~, ...
                x0_fitted, delta_x_fitted, gamma_fitted, BL_fitted, D_fitted] = ...
                fitLorentzian(bandstr(:, col_index), num_res, [], [], [], [], []);
        else
            % Initial value = parameters in the last fit
            [res_loc(:, col_index), res_wid(col_index), ~, ...
                x0_fitted, delta_x_fitted, gamma_fitted, BL_fitted, D_fitted] = ...
                fitLorentzian(bandstr(:, col_index), num_res, ...
                              x0_fitted, delta_x_fitted, gamma_fitted, BL_fitted, D_fitted);
        end
        delta_x_mean = delta_x_mean + delta_x_fitted; % Each delta_x_fitted should actually be the same
    end
    delta_x_mean = delta_x_mean / num_pts_within_roundtrip;

    real_E = zeros(1, num_pts_within_roundtrip); % unit: FSR
    real_E(1) = res_loc(round(num_res/2), 1);
    for loop_index = 2 : num_pts_within_roundtrip
        [~, temp] = min(abs(res_loc(:, loop_index) - real_E(loop_index - 1)));
        real_E(loop_index) = res_loc(temp, loop_index);
        % Find the closest resonance so that they must be on the same band.
    end
    real_E = real_E / delta_x_mean; 
    imag_E = -res_wid / delta_x_mean;
    % Renormalize the unit, from matrix index to FSR. 
    % Here real_E and imag_E are not necessarily centered to 0.
    % imag_E needs additional negative sign because of our sign convention.
    
    %% Interpolate data using periodic smoothing spline
    x = 1 : num_pts_within_roundtrip; p = 0.01; % p is smoothing parameter between 0 and 1
    f = spcsp(x, [real_E; imag_E], p);
    % This function does periodic smoothing spline.
    % https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions
    smoothed_data = fnval(f, x);
    real_E_smoothed = smoothed_data(1, :);
    imag_E_smoothed = smoothed_data(2, :);
    
    %% Finding the origin of the real and imaginary parts
    imag_E_ctr = mean(imag_E_smoothed);
    if sigma == 0, ratio = 1;
    else
        b = g_TR / (2 * kappa_TR) * cosh(sigma) / sinh(2 * sigma);
        max_f = 1 + b - 0.5;
        if b > 2, min_f = 0.5 - b; else, min_f = (-2 - b^2)/4; end
        ratio = -min_f / max_f;
    end
    real_E_ctr = min(real_E_smoothed) + ratio/(1+ratio) * range(real_E_smoothed);
    
    real_E = real_E - real_E_ctr;
    imag_E = imag_E - imag_E_ctr;
    real_E_smoothed = real_E_smoothed - real_E_ctr;
    imag_E_smoothed = imag_E_smoothed - imag_E_ctr;

    riemann_surface_real_E{sigma} = real_E;
    riemann_surface_imag_E{sigma} = imag_E;
    % Put into storage for later use
    
    %% Plot Fig. S5(B) and S5(D)
    if (sigma == 0.0) || (sigma == 0.39)
    
        figure;
        Omegat = linspace(0, 1, num_pts_within_roundtrip);
        
        yyaxis left;
        scatter(Omegat, real_E, 20, [0 0.4470 0.7410], "filled");
        ylim([-1, 1] * 0.04); yticks([-1, -0.5, 0, 0.5, 1] * 0.04);
        % ylabel('Re({\itE}) / {\it\Omega_R}');
        yticklabels([]);
        
        yyaxis right;
        scatter(Omegat, imag_E, 20, [0.8500 0.3250 0.0980], "filled"); hold on;
        ylim([-1, 1] * 0.03); yticks([-1, -0.5, 0, 0.5, 1] * 0.03);
        % ylabel('Im({\itE}) / {\it\Omega_R}');
        yticklabels([]);
        
        xlim([0, 1]); xticks([0, 0.25, 0.5, 0.75, 1]); grid off;
        % xlabel('{\itt} / {\itT_R}'); 
        xticklabels([]); 
        
        set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
            'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off');
        set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');
    end
    
    if (sigma == 0.0) || (sigma == 0.39) || (sigma == 0.5)
        %% Plot Fig. S5(F), S5(H), S5(J)
        figure; color_list = hsv(num_pts_within_roundtrip);
        scatter(real_E, imag_E, 20, color_list, 'filled');
        
        set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
            'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
            'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
        pbaspect([1, 1, 1]);
        xlim([-1, 1] * 0.04); xticks([-1, -0.5, 0, 0.5, 1] * 0.04); 
        ylim([-1, 1] * 0.03); yticks([-1, -0.5, 0, 0.5, 1] * 0.03);
        xticklabels([]); yticklabels([]); 
        % xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}');
        set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');

        %% Plot Fig. S5(G), S5(I), S5(K)
        figure; egcl = [1, 1, 1] * 0.3; fccl = [1, 1, 1] * 0.7;
        scatter(real_E, imag_E, 20, color_list, 'filled'); hold on;

        sub_reg = regions(polyshape(real_E_smoothed, imag_E_smoothed));
        switch sigma
            case 0.0
                plot(sub_reg(1), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                plot(sub_reg(2), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold off;
            case 0.39
                plot(sub_reg(1), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                plot(sub_reg(2), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                plot(sub_reg(3), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                plot(sub_reg(4), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold off;
            case 0.5
                % Winding number = 1
                plot(sub_reg(1), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                plot(sub_reg(2), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                plot(sub_reg(3), 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl); hold on;
                % Winding number = 2, we explicitly write out the coordinates
                [~, idx_11] = find_closest(sub_reg(1).Vertices(:, 1), -0.00966086);
                [~, idx_12] = find_closest(sub_reg(1).Vertices(:, 1), -0.00968435);
                [~, idx_21] = find_closest(sub_reg(2).Vertices(:, 1), -0.00968435);
                [~, idx_22] = find_closest(sub_reg(2).Vertices(:, 1), -0.0114117);
                [~, idx_31] = find_closest(sub_reg(3).Vertices(:, 1), -0.0114117);
                [~, idx_32] = find_closest(sub_reg(3).Vertices(:, 1), -0.00966086);
                sub_reg_4 = polyshape( ...
                    [sub_reg(1).Vertices(idx_12 : idx_11, 1);
                     sub_reg(3).Vertices(idx_32 : idx_31, 1); ...
                     sub_reg(2).Vertices([idx_22:end, 1:idx_21], 1); ...
                    ], ...
                    [sub_reg(1).Vertices(idx_12 : idx_11, 2); ...
                     sub_reg(3).Vertices(idx_32 : idx_31, 2); ...
                     sub_reg(2).Vertices([idx_22:end, 1:idx_21], 2); ...
                     ]);
                plot(sub_reg_4, 'LineWidth', 1, 'EdgeColor', egcl, 'FaceColor', fccl/3); hold off;
        end
        
        set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', 'TickDir', 'out', ...
            'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
            'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
        pbaspect([1, 1, 1]);
        if sigma == 0.0
            xlim([-1, 1] * 0.004); xticks([-1, 0, 1] * 0.004);
        else
            xlim([-1, 1] * 0.004 - 0.0096); xticks([-1, 0, 1] * 0.004 - 0.0096);
        end
        ylim([-1, 1] * 0.008); yticks([-1, 0, 1] * 0.008);
        xticklabels([]); yticklabels([]); 
        % xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}');
        set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.4], 'Color', 'w');
    end
end


%% Plot Riemann surface, Fig. S5(E)
figure;
for sigma = sigma_list
    color_list = hsv(num_pts_within_roundtrip);

    scatter3( ...
        riemann_surface_real_E{sigma}, ...
        riemann_surface_imag_E{sigma}, ...
        sigma * ones(1, num_pts_within_roundtrip), ...
        15, color_list, 'filled' ...
        ); hold on;

    if sigma > 0
        % Apply symmetry
        color_list = [color_list; color_list(1 : round(num_pts_within_roundtrip/2), :)];
        color_list(1 : round(num_pts_within_roundtrip/2), :) = [];
        scatter3( ...
            -riemann_surface_real_E{sigma}, ...
            -riemann_surface_imag_E{sigma}, ...
            -sigma * ones(1, num_pts_within_roundtrip), ...
            15, flipud(color_list), 'filled'); hold on;
    end

    constantplane([0, 0, 1], 0, ...
        'FaceColor', [1, 1, 1] * 0.75, 'FaceAlpha', 0.6); hold on;
    constantplane([0, 0, 1], 0.39, ...
        'FaceColor', [1, 1, 1] * 0.75, 'FaceAlpha', 0.6); hold on;

    set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
        'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'box', 'off', ...
        'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
    xlim([-1, 1] * 0.04);  xticks([-1, 0, 1] * 0.04);
    ylim([-1, 1] * 0.03);  yticks([-1, 0, 1] * 0.03);
    zlim([-1, 1] * 0.7);   zticks([-1, -0.5, 0, 0.5, 1] * 0.7);
    xticklabels([]); yticklabels([]); zticklabels([]);
    % xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}'); zlabel('Im({\itk})');
    grid off; view([-38.2, 11]);
    set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.45 0.8], 'color', 'w');
end


%% Plot OBC spectrum, Fig. S5(L)
load('./data/OBC_GBZ_simulation.mat', 'imag_k', 'real_E', 'imag_E');
% Load data; values determined manually from the winding intersections

% Apply symmetry
imag_k = [imag_k; -imag_k];
real_E = [real_E; -real_E];
imag_E = [imag_E; -imag_E];

figure; my_colormap = abyss;

% Plot simulation OBC spectrum
for plot_index = 1 : numel(real_E)
    scatter(real_E(plot_index), imag_E(plot_index), 100, ...
        my_colormap(round((imag_k(plot_index) + 0.42) / 0.84 * 255 + 1), :), 'filled'); hold on;
end
hold on;

% Plot theoretical OBC spectrum
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

real_E_rescale_factor = 0.036; imag_E_rescale_factor = 0.03;
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

pbaspect([1, 1, 1]);
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off', ...
    'YColor', [0.8500 0.3250 0.0980], 'XColor', [0 0.4470 0.7410]);
xlim([-1, 1] * 0.02); xticks([-1, 0, 1] * 0.02);
ylim([-1, 1] * 0.008); yticks([-1, 0, 1] * 0.008);
% xlabel('Re({\itE}) / {\it\Omega_R}'); ylabel('Im({\itE}) / {\it\Omega_R}');
xticklabels([]); yticklabels([]);
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.8], 'Color', 'w');

%% Plot GBZ, Fig. S5(M)
load('./data/OBC_GBZ_simulation.mat', 'imag_k', 'real_k_1', 'real_k_2');
% Load data; values determined manually from the winding intersections

% Apply symmetry
real_k = [real_k_1; real_k_2];
imag_k = [imag_k; imag_k];

z = exp(1i * (real_k + 1i * imag_k));
z = [z; -1./z];
imag_k = [imag_k; -imag_k];

figure; my_colormap = abyss;

% Plot simulation GBZ
for plot_index = 1 : numel(z)
    scatter(real(z(plot_index)), imag(z(plot_index)), 100, ...
        my_colormap(round((imag_k(plot_index) + 0.42) / 0.84 * 255 + 1), :), 'filled'); hold on;
end
hold on;

% Plot theoretical GBZ
load("./data/gbz_theory.mat"); z = Expression1;
plot(real(z), imag(z), 'LineWidth', 1, 'Color', 'k'); hold on;
line([real(z(1)), real(z(end))], [imag(z(1)), imag(z(end))], ...
    'LineWidth', 1, 'Color', 'k');

% Plot unit circle
rectangle('Position', [-1, -1, 2, 2], 'Curvature', [1, 1], ...
    'LineWidth', 1, 'LineStyle', '--');

axis equal;
set(gca, 'fontname', 'Arial', 'fontsize', 22, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', 1, 'Layer', 'Top', 'Box', 'off');
xlim([-1, 1] * 2); xticks([-2, -1, 0, 1, 2]);
ylim([-1, 1] * 2); yticks([-2, -1, 0, 1, 2]);
% xlabel('Re({\itz})'); ylabel('Im({\itz})');
xticklabels([]); yticklabels([]);
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.8], 'Color', 'w');



function [res_loc_col, res_wid, fit_lorentzian_radj, ...
    x0_fitted, delta_x_fitted, gamma_fitted, BL_fitted, D_fitted] = ...
    fitLorentzian(transmission_col, num_res_est, x0_stpt, delta_x_stpt, gamma_stpt, BL_stpt, D_stpt)
% x0: location of the first resonance peak
% delta_x: interval of the resonance peaks
% gamma: width of the resonance peaks
% BL: transmission baseline away from the resonances
% D: relative depth of the resonance peaks
expr = fittype('lorentzian_fittype(x, x0, delta_x, gamma, BL, D, num_res_loc_int)', ...
               'coefficients', {'x0', 'delta_x', 'gamma', 'BL', 'D'}, ...
               'independent', {'x'}, 'problem', {'num_res_loc_int'});
coeffName = coeffnames(expr);
stpt = zeros(1, 5); lb = zeros(1, 5);

if isempty(x0_stpt)
    stpt(strcmp('x0', coeffName)) = numel(transmission_col) * 0.125;
    stpt(strcmp('delta_x', coeffName)) = numel(transmission_col) * 0.625;
    stpt(strcmp('gamma', coeffName)) = numel(transmission_col) * 0.125;
    stpt(strcmp('BL', coeffName)) = 1;
    stpt(strcmp('D', coeffName)) = 0.5;
else
    stpt(strcmp('x0', coeffName)) = x0_stpt;
    stpt(strcmp('delta_x', coeffName)) = delta_x_stpt;
    stpt(strcmp('gamma', coeffName)) = gamma_stpt;
    stpt(strcmp('BL', coeffName)) = BL_stpt;
    stpt(strcmp('D', coeffName)) = D_stpt;
end

lb(strcmp('x0', coeffName)) = 0;
lb(strcmp('delta_x', coeffName)) = 0;
lb(strcmp('gamma', coeffName)) = 0;
lb(strcmp('BL', coeffName)) = 0;
lb(strcmp('D', coeffName)) = 0;

[f, gof] = fit(transpose(1:numel(transmission_col)), transmission_col, ...
               expr, 'problem', round(num_res_est), 'StartPoint', stpt, 'Lower', lb, 'TolFun', 1e-10);

coeffName = coeffnames(f); coeffVal = coeffvalues(f);
x0_fitted = coeffVal(strcmp('x0', coeffName));
delta_x_fitted = coeffVal(strcmp('delta_x', coeffName));
gamma_fitted = coeffVal(strcmp('gamma', coeffName));
BL_fitted = coeffVal(strcmp('BL', coeffName));
D_fitted = coeffVal(strcmp('D', coeffName));

res_loc_col = x0_fitted + (0 : (round(num_res_est) - 1))' * delta_x_fitted;
res_wid = gamma_fitted;
fit_lorentzian_radj = gof.adjrsquare;
end