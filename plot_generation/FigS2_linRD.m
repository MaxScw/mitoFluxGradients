%% load setup file for file paths and data
clc
clear all
temp_resolved = false;
run("setup.m")

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
xlabel('distance from cell center ($r$) [$\mu$  m]', 'Interpreter','latex');
ylabel('$J_{\mathrm{ox}}$ [$\mu$  M/s]', 'Interpreter','latex');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',19);

% cb = colorbar;
% clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
% title(cb, 'c_{oxy} ($\mu$  M)', 'Interpreter', 'tex')
% cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);
xlim([4, 36])
ylim([-25, 140])

cb = colorbar;
clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
title(cb, '$c_{\mathrm{out}}$ ($\mu$  M)', 'Interpreter', 'latex')
cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);

savefig(string(plot_path)+'Jox_vs_dist_fit_linRD'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_fit_linRD'+string(temp_string)+'.png')
%%
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');
errorbar(precise_oxygen_levels, fit_params(2, :), sigma_fit_params(2, :), 'o', 'LineWidth',1.5, 'MarkerSize',10)
set(gca,'FontSize',19);
xlim([-5 55])
xticks([0 10 20 30 40 50])
xlabel("$c_\mathrm{out}$ [$\mu$  M]", 'Interpreter','latex')
ylabel("amplitude $A$ [$\mu$  M/s]", 'Interpreter','latex')
box off

savefig(string(plot_path)+'Jox_vs_dist_fit_linRD_paramAmp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_fit_linRD_paramAmp'+string(temp_string)+'.png')
%%
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');
errorbar(precise_oxygen_levels, fit_params(3, :), sigma_fit_params(3, :), 'o', 'LineWidth',1.5, 'MarkerSize',10)
set(gca,'FontSize',19);
xlim([-5 55])
xticks([0 10 20 30 40 50])
xlabel("$c_\mathrm{out}$ [$\mu$  M]", 'Interpreter','latex')
ylabel("decay length $\lambda$ [$\mu$  m]", 'Interpreter','latex')
box off

savefig(string(plot_path)+'Jox_vs_dist_linRD_paramDec.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_linRD_paramDec.png')