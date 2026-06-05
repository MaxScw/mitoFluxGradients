%% designation of oxygen ranges to analyse
clear all;
clc;

oxygen_ranges={'run_oxy_l0_h0p1', 'run_oxy_l0p1_h0p2', 'run_oxy_l0p2_h0p3', 'run_oxy_l0p3_h0p4',...
               'run_oxy_l0p4_h0p6', 'run_oxy_l0p6_h0p8',...
               'run_oxy_l0p8_h1', 'run_oxy_l1_h1p2',...
               'run_oxy_l1p2_h1p4', 'run_oxy_l1p4_h1p6',...
               'run_oxy_l1p6_h1p8',...
               'run_oxy_l1p8_h2', 'run_oxy_l2_h2p8',...
               'run_oxy_l2p8_h3p6', 'run_oxy_l3p6_h4p4',...
               'run_oxy_l4p4_h5p2'};
num_oxygen_ranges={[0, 0.1], ...
                   [0.1, 0.2], [0.2, 0.3], [0.3, 0.4],...
                   [0.4, 0.6], [0.6, 0.8], [0.8, 1.0], [1.0, 1.2],...
                   [1.2, 1.4], [1.4, 1.6], [1.6, 1.8], [1.8, 2.0],...
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4], [4.4, 5.2]};

temp_resolved = false;
data_path = "/home/mx/mitoFluxGradients/data/";
pp_path = '../data/EXP_published/postprocessing_results/';
plot_path = '../data/EXP_published/plots/';
% temp_resolved = true;
% data_path = "/home/mx/mitoFluxGradients/data/";
% pp_path = '../data/EXP_temperatureResolved/postprocessing_results/';
% plot_path = '../data/EXP_temperatureResolved/plots/';

load(data_path + 'exp_solubilities.mat')
load(data_path + 'exp_temperatures.mat')

%% loading of mito density, temperatures, solubilities, diffusivities
load(data_path + "mito_density_oocytes/mito_density.mat")
mean_mito_density = mean(ratio_num_all_mitotracker, 1);

load(data_path + "exp_solubilities.mat")
load(data_path + "exp_temperatures.mat")
load(data_path + "exp_diffusivities.mat")

if temp_resolved==true
    temp_ind = 4;
    temp_string = '_T'+string(temperature(temp_ind))+'C';
else 
    temp_string = '';
end

%% loading of read-in (experimental) data

% choose which temperature to plot analysis for (from T=[22, 28, 31, 36]C)
load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');

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
sigma_precise_oxygen_levels = (213.5/20.946).*sigma_precise_oxygen_levels;
sigma_precise_oxygen_levels(sigma_precise_oxygen_levels==0) = mean(sigma_precise_oxygen_levels);
end


%%
load(pp_path + "gradient_fit_linear_model"+string(temp_string)+".mat")
fit_params = fit_data_gradient{1};
sigma_fit_params = fit_data_gradient{2};

%%

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],'color','white');

cs = flip(viridis(20));
colormap(cs)

start_ring = 2;


for ind=1:16
    
    hold on
    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    % errorbar(data(1, :), data(2, :).*mean_mito_density, data_stderr(2, :), data_stderr(2, :),...
    %          'o', 'Color', cs(ind, :))
    
 

    % USE MITO-DENSITY WEIGHTED JOX DATA
    mito_weighted_data = data(2, :);%0.5*data(2, :).*mean(mean_mito_density);
    
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    mito_weighted_data_stderr = 0.5*mito_weighted_data.*sqrt((data_stderr(2, :)./data(2, :)).^2);
    if ind==16
        errorbar(data(1, start_ring:end), mito_weighted_data(start_ring:end), mito_weighted_data_stderr(start_ring:end), mito_weighted_data_stderr(start_ring:end),...
                 'o', 'Color', cs(ind, :), 'LineWidth',1.5, 'DisplayName','$J_\mathrm{ox}$')
        plot(data(1, :), rd_solver_linear(data(1, :), fit_params(:, ind)), 'Color',cs(ind, :), ...
            'LineWidth',1.5, 'DisplayName', 'lin. RD fit')
    else
       errorbar(data(1, start_ring:end), mito_weighted_data(start_ring:end), mito_weighted_data_stderr(start_ring:end), mito_weighted_data_stderr(start_ring:end),...
                 'o', 'Color', cs(ind, :), 'LineWidth',1.5, 'HandleVisibility', 'off')
        plot(data(1, start_ring:end), rd_solver_linear(data(1, start_ring:end), fit_params(:, ind)), ...
            'Color',cs(ind, :), 'LineWidth',1.5, 'HandleVisibility', 'off')
    end
end

legend(Interpreter="latex")
xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex');
ylabel('$J_{\mathrm{ox}}$ [\,\,\,\,M/s]', 'Interpreter','latex');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',19);

% cb = colorbar;
% clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
% title(cb, 'c_{oxy} (\mu M)', 'Interpreter', 'tex')
% cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);
xlim([4, 36])
ylim([-25, 140])

cb = colorbar;
clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
title(cb, '$c_{\mathrm{out}}$ (\,\,\,\,M)', 'Interpreter', 'latex')
cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);

savefig(string(plot_path)+'Jox_vs_dist_fit_linRD'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_fit_linRD'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'Jox_vs_dist_fit_linRD'+string(temp_string)+'.eps')
%%
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');
errorbar(precise_oxygen_levels, fit_params(2, :), sigma_fit_params(2, :), 'o', 'LineWidth',1.5, 'MarkerSize',10)
set(gca,'FontSize',19);
xlim([-5 55])
xticks([0 10 20 30 40 50])
xlabel("$c_\mathrm{out}$ [\,\,\,\,M]", 'Interpreter','latex')
ylabel("amplitude $A$ [\,\,\,\,M/s]", 'Interpreter','latex')
box off

savefig(string(plot_path)+'Jox_vs_dist_fit_linRD_paramAmp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_fit_linRD_paramAmp'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'Jox_vs_dist_fit_linRD_paramAmp'+string(temp_string)+'.eps')
%%
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');
errorbar(precise_oxygen_levels, fit_params(3, :), sigma_fit_params(3, :), 'o', 'LineWidth',1.5, 'MarkerSize',10)
set(gca,'FontSize',19);
xlim([-5 55])
xticks([0 10 20 30 40 50])
xlabel("$c_\mathrm{out}$ [\,\,\,\,M]", 'Interpreter','latex')
ylabel("decay length $\lambda$ [\,\,\,\,m]", 'Interpreter','latex')
box off

savefig(string(plot_path)+'Jox_vs_dist_linRD_paramDec.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_linRD_paramDec.png')
export_fig(fig, string(plot_path)+'Jox_vs_dist_linRD_paramDec'+string(temp_string)+'.eps')