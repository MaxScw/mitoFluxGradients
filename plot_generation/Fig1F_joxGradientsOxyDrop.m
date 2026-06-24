%% load setup file for file paths and data
clc
clear all
temp_resolved = false;
run("setup.m")

%% simple plot of J_ox gradients for different outside oxygen

% subplot with all oxygen levels included
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');

cs = flip(viridis(20));
colormap(cs)


%title(temp_string, 'Interpreter','none')

hold on;
for k=1:numel(oxygen_ranges)
    
    dist_all=oxygen_ranges_data{k}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{k}.jox_cell_kn_all;

    errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(jox_cell_kn_all(:, 2:end), 'omitnan'),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));

    xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex');
    ylabel('$J_{\mathrm{ox}}$ [$\mu$ M/s]', 'Interpreter','latex');
    
    ylim([-18, 140])

    set(gca,'FontSize',19);
end

%ylim([-10, 140])
xlim([4, 36])


cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$c_\mathrm{out}$ [$\mu$ M]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);




% subplot with 3 maximally different levels

axes(Position=[.29 .7 .2 .2])
box on

%fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 400]);
cs = flip(viridis(20));
colormap(cs)

%title(temp_string, 'Interpreter','none')
hold on
for k=[1, 2, numel(oxygen_ranges)]
    
    dist_all=oxygen_ranges_data{k}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{k}.jox_cell_kn_all;

    errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(jox_cell_kn_all(:, 2:end), 'omitnan'),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));

    %xlabel('distance from cell center $(\mu m)$', 'Interpreter','latex');
    %ylabel('$J_{ox} (\mu M/s)$', 'Interpreter','latex');
    
    ylim([-18, 100])

    set(gca,'FontSize',19);
end

%ylim([-10, 140])
xlim([4, 36])
%ylim([-25, 140])


savefig(string(plot_path)+'Jox_vs_dist_emp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_emp'+string(temp_string)+'.png')