%% designation of oxygen ranges to analyse
clear all;
clc;

oxygen_ranges={"run_oxy_l0_h0p1", "run_oxy_l0p1_h0p2", "run_oxy_l0p2_h0p3", "run_oxy_l0p3_h0p4",...
               "run_oxy_l0p4_h0p6", "run_oxy_l0p6_h0p8",...
               "run_oxy_l0p8_h1", "run_oxy_l1_h1p2",...
               "run_oxy_l1p2_h1p4", "run_oxy_l1p4_h1p6",...
               "run_oxy_l1p6_h1p8",...
               "run_oxy_l1p8_h2", "run_oxy_l2_h2p8",...
               "run_oxy_l2p8_h3p6", "run_oxy_l3p6_h4p4",...
               "run_oxy_l4p4_h5p2"};
num_oxygen_ranges={[0, 0.1], ...
                   [0.1, 0.2], [0.2, 0.3], [0.3, 0.4],...
                   [0.4, 0.6], [0.6, 0.8], [0.8, 1.0], [1.0, 1.2],...
                   [1.2, 1.4], [1.4, 1.6], [1.6, 1.8], [1.8, 2.0],...
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4],... 
                   [4.4, 5.2]};

temp_resolved = false;
data_path = "/home/mx/mitoFluxGradients/data/";
pp_path = "/home/mx/mitoFluxGradients/data/EXP_published/postprocessing_results/";
plot_path = "/home/mx/mitoFluxGradients/data/EXP_published/plots/";
% temp_resolved = true;
% data_path = "/home/mx/mitoFluxGradients/data/";
% pp_path = "../data/EXP_temperatureResolved/postprocessing_results/";
% plot_path = "../data/EXP_temperatureResolved/plots/";

%% loading of mito density, temperatures, solubilities, diffusivities

load(data_path + "/mito_density_oocytes/mito_density.mat")
mean_mito_density = mean(ratio_num_all_mitotracker, 1);
sigma_mito_weigth = std(ratio_num_all_mitotracker, 1)./numel(ratio_num_all_mitotracker(1, :));

load(data_path + "exp_solubilities.mat")
load(data_path + "exp_temperatures.mat")
load(data_path + "exp_diffusivities.mat")

if temp_resolved==true
    temp_ind = 1;
    temp_string = "_T"+string(temperature(temp_ind))+"C"
else
    temp_ind = 4;
    temp_string = "";
end

%% loading of read-in data

load(pp_path + "plot_data_multiple_oxy_ranges" + temp_string + ".mat");

precise_oxygen_levels = [];
sigma_precise_oxygen_levels = [];
for oxy_ind=1:numel(oxygen_ranges)
    
    precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
    sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, mean(oxygen_ranges_data{oxy_ind}.sigma_o2_levels, 'omitnan')];
end


if temp_resolved==true
precise_oxygen_levels = solubility(temp_ind).*precise_oxygen_levels./20.946;
sigma_precise_oxygen_levels = solubility(temp_ind).*sigma_precise_oxygen_levels./20.946;
sigma_precise_oxygen_levels(sigma_precise_oxygen_levels==0) = mean(sigma_precise_oxygen_levels);
else
precise_oxygen_levels = (213.5/20.946).*precise_oxygen_levels;
end

%% fit gradient data to different reaction-diffusion models

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];

% get mean data + do fits
for ind=6%:numel(oxygen_ranges)
    
    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
        data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];

    % USE SPACE-AVERAGED OR INHOM. MITO WEIGHT
    mito_weight = mean(mean_mito_density);

    % MITO-DENSITY WEIGHTED JOX DATA
    mito_weighted_data = data(2, :).*mito_weight.*0.5;
    mito_weighted_data_stderr = 0.5*mito_weighted_data.*sqrt((data_stderr(2, :)./data(2, :)).^2 + (sigma_mito_weigth./mito_weight).^2);
    
    %use mito weighted data
    data(2, :) = mito_weighted_data;
    data_stderr(2, :) = mito_weighted_data_stderr;
    
    y_weights = 1./data_stderr(2, :);
    
    
    % set parameters to fit, increment and initial values (MICHAELIS-MENTEN RD, fixed km)
    % param = [c_star, v_max, k_m, D]
    
    dp = [0 0.001 0.001 0];

    km0 = data(2, end);
    vm0 = data(2, end)*100;
    p0 = [precise_oxygen_levels(ind) vm0 km0 diffusivity(temp_ind)];

    pmin = [0.1 0.0001 0.0001 2000];
    pmax = [200 1e6 1e6 3000];

    % fit
    
    fit_start = 2;
    fit_end = numel(data(1, :))
    max_iter = 1000;
    
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@rd_solver_mm3d, p0,...
                     data(1, :)', data(2, :)', y_weights(fit_start:fit_end), dp,...
                     pmin, pmax, [], fit_start, fit_end, max_iter);
    
    % !!!optional:
    % use sum of squared residuals instead of built-in method of lm to evaluate X2
    resid = (data(2, fit_start:fit_end) - rd_solver_mm3d(data(1, fit_start:fit_end), p)');
    
    x2_jox = sum((resid).^2./(mito_weighted_data(fit_start:fit_end)));
    X2 = x2_jox;
    
    R_sq = corrcoef(data(2, fit_start:fit_end), rd_solver_mm3d(data(1, fit_start:fit_end), p)');
    R_sq = R_sq(1,2).^2;

    fit_params = [fit_params, p];
    sigma_fit_params = [sigma_fit_params, sigma_p];
    chi_squareds = [chi_squareds X2];
    R_squareds = [R_squareds R_sq];

    plot(data(1, :), data(2, :), 'o')
    hold on
    plot(data(1, :), rd_solver_mm3d(data(1, :), p0), 'Color','red')
    plot(data(1, :), rd_solver_mm3d(data(1, :), p), 'Color','green')
end

fit_data_gradient{1} = fit_params;
fit_data_gradient{2} = sigma_fit_params;
fit_data_gradient{3} = chi_squareds;
fit_data_gradient{4} = R_squareds;

save(pp_path + "gradient_fit_mm_num_sol_model"+string(temp_string)+".mat", 'fit_data_gradient');
%%
mean(chi_squareds)
mean(R_squareds)
%%
subplot(1, 3, 1)
errorbar(fit_params(2, :), sigma_fit_params(2, :), sigma_fit_params(2, :), 'o')
subplot(1, 3, 2)
errorbar(fit_params(3, :), sigma_fit_params(3, :), sigma_fit_params(3, :), 'o')
subplot(1, 3, 3)
errorbar(fit_params(4, :), sigma_fit_params(4, :), sigma_fit_params(4, :), 'o')

figure(2)
plot(R_squareds, 'o', 'HandleVisibility','off')
yline(mean(R_squareds), 'DisplayName','\langle R^2 \rangle='+string(mean(R_squareds)))
legend()

%%

figure('Renderer', 'painters', 'Position', [10 10 1500 500])

cs = viridis(16);
colormap(cs)
hold on

for ind=1:6
    subplot(1, 3, 1)
    hold on
    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    c = data(2, end);

    % USE SPACE-AVERAGED OR INHOM. MITO WEIGHT
    mito_weight = mean(mito_weight);

    % MITO-DENSITY WEIGHTED JOX DATA
    mito_weighted_data = data(2, :).*mito_weight.*0.5;
    mito_weighted_data_stderr = 0.5*mito_weighted_data.*sqrt((data_stderr(2, :)./data(2, :)).^2 + (sigma_mito_weigth./mito_weight).^2);
    
    %use mito weighted data
    data(2, :) = mito_weighted_data;
    data_stderr(2, :) = mito_weighted_data_stderr;

    % errorbar(data(1, :), data(2, :).*mito_weight, data_stderr(2, :), data_stderr(2, :),...
    %          'o', 'Color', cs(ind, :))
    errorbar(data(1, :), mito_weighted_data, data_stderr(2, :), data_stderr(2, :),...
             'o', 'Color', cs(ind, :), 'LineWidth',1.5)
    plot(data(1, :), rd_solver_mm3d(data(1, :), fit_params(:, ind), c), 'Color',cs(ind, :), 'LineWidth',1.5)
end
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',15);

% cb = colorbar;
% clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
% title(cb, 'c_{oxy} (\mu M)', 'Interpreter', 'tex')
% cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);

for ind=7:12
    subplot(1, 3, 2)
    hold on
    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    c = data(2, end);

    % USE SPACE-AVERAGED OR INHOM. MITO WEIGHT
    mito_weight = mean(mito_weight);

    % MITO-DENSITY WEIGHTED JOX DATA
    mito_weighted_data = data(2, :).*mito_weight.*0.5;
    mito_weighted_data_stderr = 0.5*mito_weighted_data.*sqrt((data_stderr(2, :)./data(2, :)).^2 + (sigma_mito_weigth./mito_weight).^2);
    
    %use mito weighted data
    data(2, :) = mito_weighted_data;
    data_stderr(2, :) = mito_weighted_data_stderr;

    % errorbar(data(1, :), data(2, :).*mito_weight, data_stderr(2, :), data_stderr(2, :),...
    %          'o', 'Color', cs(ind, :))
    errorbar(data(1, :), data(2, :), data_stderr(2, :), data_stderr(2, :),...
             'o', 'Color', cs(ind, :), 'LineWidth',1.5)
    plot(data(1, :), rd_solver_mm3d(data(1, :), fit_params(:, ind), c), 'Color',cs(ind, :), 'LineWidth',1.5)
end
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',15);
% 
% cb = colorbar;
% clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
% title(cb, 'c_{oxy} (\mu M)', 'Interpreter', 'tex')
% cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);

for ind=13:16
    subplot(1, 3, 3)
    hold on
    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    c = data(2, end);

    % USE SPACE-AVERAGED OR INHOM. MITO WEIGHT
    mito_weight = mean(mito_weight);

    % MITO-DENSITY WEIGHTED JOX DATA
    mito_weighted_data = data(2, :).*mito_weight.*0.5;
    mito_weighted_data_stderr = 0.5*mito_weighted_data.*sqrt((data_stderr(2, :)./data(2, :)).^2 + (sigma_mito_weigth./mito_weight).^2);
    
    %use mito weighted data
    data(2, :) = mito_weighted_data;
    data_stderr(2, :) = mito_weighted_data_stderr;

    % errorbar(data(1, :), data(2, :).*mito_weight, data_stderr(2, :), data_stderr(2, :),...
    %          'o', 'Color', cs(ind, :))
    errorbar(data(1, :), data(2, :), data_stderr(2, :), data_stderr(2, :),...
             'o', 'Color', cs(ind, :), 'LineWidth',1.5)
    plot(data(1, :), rd_solver_mm3d(data(1, :), fit_params(:, ind), c), 'Color',cs(ind, :), 'LineWidth',1.5)
end
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',15);

cb = colorbar;
clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
title(cb, 'c_{oxy} (\mu M)', 'Interpreter', 'tex')
cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);
