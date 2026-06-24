%% load setup file for file paths and data
clc
clear all
temp_resolved = false;
run("setup.m")



%% plot Jox as a function of subcellular oxygen for all rings

jox_per_ring = [];
sigma_jox_per_ring = [];
dist_per_ring = [];

for ind=1:numel(oxygen_ranges)

    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    jox_per_ring = [jox_per_ring, data(2, :)'];
    sigma_jox_per_ring = [sigma_jox_per_ring, data_stderr(2, :)'];
    dist_per_ring = [dist_per_ring, data(1, :)'];
end

% choose all or only selection of all distances (rings) to plot
plot_selectedDist = false;
dists = [2, 8, 10];

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600], 'color','white');
% save corrected oxygen levels
load(string(pp_path)+'cOxy_corr'+string(temp_string)+'.mat')

fs = 18;

% subplot(1, 2, 2)
hold on;
colors = flip(cool(15));
colormap(colors(1:10, 1:end))

xlabel('$\hat{c}(c_\mathrm{out}, r)$ [$\mu$ M]', 'Interpreter','latex');
ylabel('$J_\mathrm{ox}$ [$\mu$ M/s]', 'Interpreter','latex');
xlim([-5 55])
ylim([-4 135])

set(gca,'FontSize',fs);

cb = colorbar;
clim([0.5 10.5])
title(cb, '$r$ [$\mu$ m]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, 10, 10);
cb.TickLabels = round(mean(dist_per_ring, 2), 2);
title(temp_string, 'Interpreter','none')

if plot_selectedDist==false
    for k=2:10
        for j=1:numel(oxygen_ranges)
            dist_all=oxygen_ranges_data{j}.dist_all;
            jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
            if (j==1)&&(k==10)
                errorbar(cOxy_all(2, j, k),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                     std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                     'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',colors(k, :), ...
                     'DisplayName','$J_\mathrm{ox}$');
            else
                errorbar(cOxy_all(2, j, k),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                         std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                         'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',colors(k, :), ...
                         'HandleVisibility','off');
            end
            
        end
    end

    
    % axes('Position',[.6 .7 .2 .2])
    % hold on
    % yscale("log")
    % %ylabel("log($J_\mathrm{ox}$)", Interpreter='latex', FontSize=16)
    % for k=2:10
    %     for j=1:4
    %         dist_all=oxygen_ranges_data{j}.dist_all;
    %         jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
    %         errorbar(squeeze(cOxy_all(2, j, k)),mean(jox_cell_kn_all(:, k), 'omitnan'),...
    %              std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1))./ ...
    %              mean(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
    %              'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',colors(k, :));
    %     end    
    % end
else
    for k=dists
        for j=1:numel(oxygen_ranges)
            dist_all=oxygen_ranges_data{j}.dist_all;
            jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
            errorbar(cOxy_all(2, j, k),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                     std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                     'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',colors(k, :));
        end
    end

    % box on
    % axes('Position',[.6 .7 .2 .2])
    % hold on
    % for k=dists
    %     for j=1:7
    %         dist_all=oxygen_ranges_data{j}.dist_all;
    %         jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
    %         errorbar(squeeze(cOxy_all(2, j, k)),mean(jox_cell_kn_all(:, k), 'omitnan'),...
    %              std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
    %              'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',colors(k, :));
    %     end    
    %end
end
%legend('Interpreter','latex')


if plot_selectedDist==false
    savefig(string(plot_path)+'Jox_vs_coxy_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'Jox_vs_coxy_emp'+string(temp_string)+'.png')
else
    savefig(string(plot_path)+'Jox_vs_coxy_selected_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'Jox_vs_coxy_selected_emp'+string(temp_string)+'.png')
end
