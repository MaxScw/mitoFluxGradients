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

%% load in single-temp data
temp_ind = 4;
    
if temp_resolved==true
    temp_string = '_T'+string(temperature(temp_ind))+'C'
else
    temp_ind = 1;
    temp_string = '';
end
load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');

precise_oxygen_levels = [];
sigma_precise_oxygen_levels = [];
for oxy_ind=1:numel(oxygen_ranges)
    precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
    sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.sigma_o2_levels];
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


%% plot highest oxy level jox gradient for all temps
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],...
    'Color','white');
cs = plasma(4+1);
cs = cs(1:end-1, 1:end);
colormap(cs)

hold on

xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex');
ylabel('$J_{\textrm{ox}}$ [\,\,\,\,M/s]', 'Interpreter','latex');

set(gca,'FontSize',19);

cb = colorbar;
clim([0.5 4.5])
title(cb, 'T [$^\circ$C]', 'Interpreter', 'latex')

cb.Ticks = linspace(1, 4, 4);
kelvin_temps = temperature(1:end-1)+273.15;
cb.TickLabels = temperature(1:end-1);

max_oxy = 50; % muM
min_oxy = 20;
title_string = 'average over '+string(min_oxy)+'$< c_\mathrm{out} <$'+string(max_oxy)+'\,\,\,\,M';
%title(title_string, 'Interpreter','latex')

xlim([0, 36])

for t=1:4

    temp_ind = t;
    
    if temp_resolved==true
        temp_string = '_T'+string(temperature(temp_ind))+'C';
    else
        temp_ind = 4;
        temp_string = '';
    end
    load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
    
    precise_oxygen_levels = [];
    sigma_precise_oxygen_levels = [];
    for oxy_ind=1:numel(oxygen_ranges)
        precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
        sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.sigma_o2_levels];
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

    % plotting
    
    dist_av = [];
    jox_av = [];
    jox_std_av = [];

    for oxy=1:16
        if precise_oxygen_levels(oxy)<=max_oxy
            if precise_oxygen_levels(oxy)>=min_oxy
                dist_all=oxygen_ranges_data{oxy}.dist_all;
                jox_cell_kn_all=oxygen_ranges_data{oxy}.jox_cell_kn_all;

                dist_av = [dist_av; mean(dist_all(:, 2:end), 'omitnan')];
                jox_av = [jox_av; mean(jox_cell_kn_all(:, 2:end), 'omitnan')];
                jox_std_av = [jox_std_av; std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1))];
            end
        end
    end
    
    jox = mean(jox_av, 1);
    jox_std = mean(jox_std_av, 1);
    dist = mean(dist_av, 1);

    dist_all=oxygen_ranges_data{end}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{end}.jox_cell_kn_all;

    % errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(jox_cell_kn_all(:, 2:end), 'omitnan'),...
    %          std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
    %          std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
    %          'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(t, :));
    errorbar(dist, jox, jox_std, jox_std, 'o', 'MarkerSize',15,'LineWidth',1.5, 'Color',cs(t, :));

    %ylim([-20, 80])
    

end

% subplot with 3 maximally different levels

axes(Position=[.31 .75 .2 .2])
box on
hold on

ylim([0.35 1.1])
xlim([0 37])
xticks([0 10 20 30])
yticks([0.5 0.75 1])
set(gca,'FontSize',19);

for t=1:4

    temp_ind = t;
    
    if temp_resolved==true
        temp_string = '_T'+string(temperature(temp_ind))+'C';
    else
        temp_ind = 4;
        temp_string = '';
    end
    load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
    
    precise_oxygen_levels = [];
    sigma_precise_oxygen_levels = [];
    for oxy_ind=1:numel(oxygen_ranges)
        precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
        sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.sigma_o2_levels];
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

    % plotting
    
    dist_av = [];
    jox_av = [];
    jox_std_av = [];

    for oxy=1:16
        if precise_oxygen_levels(oxy)<=max_oxy
            if precise_oxygen_levels(oxy)>=min_oxy
                dist_all=oxygen_ranges_data{oxy}.dist_all;
                jox_cell_kn_all=oxygen_ranges_data{oxy}.jox_cell_kn_all;

                dist_av = [dist_av; mean(dist_all(:, 2:end), 'omitnan')];
                jox_av = [jox_av; mean(jox_cell_kn_all(:, 2:end), 'omitnan')];
                jox_std_av = [jox_std_av; std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1))];
            end
        end
    end
    
    jox = mean(jox_av, 1);
    jox_std = mean(jox_std_av, 1);
    dist = mean(dist_av, 1);

    dist_all=oxygen_ranges_data{end}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{end}.jox_cell_kn_all;

    % errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(jox_cell_kn_all(:, 2:end), 'omitnan'),...
    %          std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
    %          std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
    %          'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(t, :));
    %errorbar(dist, jox./jox(end), jox_std./jox_std(end), jox_std./jox_std(end), 'o', 'MarkerSize',15,'LineWidth',1.5, 'Color',cs(t, :));
    plot(dist, jox./jox(end),'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(t, :));

    %ylim([-20, 80])
    

end

savefig(string(plot_path)+'Jox_vs_dist_emp_tempComp.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_emp_tempComp.png')
export_fig(fig, string(plot_path)+'Jox_vs_dist_emp_tempComp.eps')
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

    xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex');
    ylabel('$J_{\mathrm{ox}}$ [\,\,\,\,M/s]', 'Interpreter','latex');
    
    ylim([-18, 140])

    set(gca,'FontSize',19);
end

%ylim([-10, 140])
xlim([4, 36])


cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$c_\mathrm{out}$ [\,\,\,\,M]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);




% subplot with 3 maximally different levels

axes(Position=[.29 .7 .2 .2])
box on

%fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 400]);
cs = flip(viridis(20));
colormap(cs)

title(temp_string, 'Interpreter','none')
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
export_fig(fig, string(plot_path)+'Jox_vs_dist_emp'+string(temp_string)+'.eps')
% cb = colorbar;
% clim([1.5 numel(precise_oxygen_levels)+0.5])
% title(cb, '$c_\mathrm{out} (\mu M)$', 'Interpreter', 'latex')
% cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
% cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);

%savefig(string(plot_path)+'Jox_vs_dist_emp_maxDiff'+string(temp_string)+'.fig')
%saveas(fig2, string(plot_path)+'Jox_vs_dist_emp_maxDiff'+string(temp_string)+'.png')
%saveas(fig2, string(plot_path)+'Jox_vs_dist_emp_maxDiff'+string(temp_string)+'.svg')

%% plot NADH free and bound for different external oxygen for all rings

% subplot with all oxygen levels included
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');
cs = flip(viridis(20));
colormap(cs)


%title(temp_string, 'Interpreter','none')

hold on;
for k=1:numel(oxygen_ranges)
    
    dist_all=oxygen_ranges_data{k}.dist_all;
    nadhf_all=oxygen_ranges_data{k}.nadhf;
    nadhb_all=oxygen_ranges_data{k}.nadhb;

    errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(nadhf_all(:, 2:end), 'omitnan'),...
             std(nadhf_all(:, 2:end), 'omitnan')./sqrt(size(nadhf_all(:, 2:end),1)),...
             std(nadhf_all(:, 2:end), 'omitnan')./sqrt(size(nadhf_all(:, 2:end),1)),...
             'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));

    % errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(nadhb_all(:, 2:end), 'omitnan'),...
    %          std(nadhb_all(:, 2:end), 'omitnan')./sqrt(size(nadhb_all(:, 2:end),1)),...
    %          std(nadhb_all(:, 2:end), 'omitnan')./sqrt(size(nadhb_all(:, 2:end),1)),...
    %          'o', 'MarkerSize',10,'LineWidth',1.5, 'Color','red');
    % 
    % errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(nadhf_all(:, 2:end), 'omitnan'),...
    %          std(nadhf_all(:, 2:end), 'omitnan')./sqrt(size(nadhf_all(:, 2:end),1)),...
    %          std(nadhf_all(:, 2:end), 'omitnan')./sqrt(size(nadhf_all(:, 2:end),1)),...
    %          'o', 'MarkerSize',10,'LineWidth',1.5, 'Color','blue');

    xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex');
    ylabel('$[\mathrm{NADH}_\mathrm{f}]$ [\,\,\,\,M]', 'Interpreter','latex');

    

    set(gca,'FontSize',19);
end

if (temp_resolved==true)&&(temp_ind==1)
cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$c_\mathrm{out}$ [\,\,\,\,M]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);
end

title('T='+string(temperature(temp_ind))+'$^\circ C$', 'Interpreter', 'latex')

xlim([4, 36])
%ylim([20 85])


savefig(string(plot_path)+'NADHf_vs_dist_emp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'NADHf_vs_dist_emp'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'NADHf_vs_dist_emp'+string(temp_string)+'.eps')



% subplot with all oxygen levels included
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],'color','white');
cs = flip(viridis(20));
colormap(cs)


%title(temp_string, 'Interpreter','none')

hold on;
for k=1:numel(oxygen_ranges)
    
    dist_all=oxygen_ranges_data{k}.dist_all;
    nadhf_all=oxygen_ranges_data{k}.nadhf;
    nadhb_all=oxygen_ranges_data{k}.nadhb;

    
    errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(nadhb_all(:, 2:end), 'omitnan'),...
             std(nadhb_all(:, 2:end), 'omitnan')./sqrt(size(nadhb_all(:, 2:end),1)),...
             std(nadhb_all(:, 2:end), 'omitnan')./sqrt(size(nadhb_all(:, 2:end),1)),...
             'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));

    % errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(nadhf_all(:, 2:end), 'omitnan'),...
    %          std(nadhf_all(:, 2:end), 'omitnan')./sqrt(size(nadhf_all(:, 2:end),1)),...
    %          std(nadhf_all(:, 2:end), 'omitnan')./sqrt(size(nadhf_all(:, 2:end),1)),...
    %          'o', 'MarkerSize',10,'LineWidth',1.5, 'Color','blue');

    xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex');
    ylabel('$[\mathrm{NADH}_\mathrm{b}]$ [\,\,\,\,M]', 'Interpreter','latex');

    

    set(gca,'FontSize',19);
end

if (temp_resolved==true)&&(temp_ind==1)
cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$\mathrm{c}_{\mathrm{out}}$ [\,\,\,\,M]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);
end

title('T='+string(temperature(temp_ind))+'$^\circ C$', 'Interpreter', 'latex')

xlim([4, 36])
%ylim([15 40])

savefig(string(plot_path)+'NADHb_vs_dist_emp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'NADHb_vs_dist_emp'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'NADHb_vs_dist_emp'+string(temp_string)+'.eps')

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

fs = 30;

% subplot(1, 2, 2)
hold on;
colors = flip(cool(15));
colormap(colors(1:10, 1:end))

xlabel('$\hat{c}(c_\mathrm{out}, r)$ [\,\,\,\,M]', 'Interpreter','latex');
ylabel('$J_\mathrm{ox}$ [\,\,\,\,M/s]', 'Interpreter','latex');
xlim([-5 55])
ylim([-4 135])

set(gca,'FontSize',fs);

cb = colorbar;
clim([0.5 10.5])
title(cb, '$r$ [\,\,\,\,m]', 'Interpreter', 'latex')
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
    export_fig(fig, string(plot_path)+'Jox_vs_coxy_emp'+string(temp_string)+'.eps')
else
    savefig(string(plot_path)+'Jox_vs_coxy_selected_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'Jox_vs_coxy_selected_emp'+string(temp_string)+'.png')
end
%% calculate mean squared distance (+ ringw-wise covariance, kolmogorow-smirnov test) between oxy levels 

covar = zeros(numel(oxygen_ranges)-1, 10);
kstest2_results = zeros(numel(oxygen_ranges)-1, 10);
kstest2_p = zeros(numel(oxygen_ranges)-1, 10);
klDiv = zeros(numel(oxygen_ranges)-1, 10);
msds = zeros(numel(oxygen_ranges), 0);
corr_matrices = {};
p_vals = {};

for k=1:numel(oxygen_ranges)-1
    k
    current = oxygen_ranges_data{k}.jox_cell_kn_all;
    next = oxygen_ranges_data{k+1}.jox_cell_kn_all;
    trunc_size = min(size(current, 1), size(next, 1));

    msds(k) = mean((mean(current, 'omitnan') - mean(next, 'omitnan')).^2);
    
    [rho,pval] = corr(current(1:trunc_size, :), next(1:trunc_size, :), 'Type', 'Kendall', 'rows', 'complete');
    corr_matrices{k} = rho;
    p_vals{k} = pval;

    for j=1:min(size(current, 2), size(next, 2))
        [kstest2_results(k, j), kstest2_p(k, j)] = kstest2(current(:, j), next(:, j));

        dev_current = std(current(:, j), 'omitnan')./sqrt(size(current(:, j),1));
        dev_next = std(next(:, j), 'omitnan')./sqrt(size(next(:, j),1));

        cov_mat = cov(current(1:trunc_size, j), next(1:trunc_size, j), 'omitrows');
        covar(k, j) = (cov_mat(1, 2) + cov_mat(2, 1))/2;
        covar(k, j) = covar(k, j)/(dev_current*dev_next);
        
    end
end

%% create list of correlation coefficients corresponding to oxygen levels

mean_rho = zeros(numel(oxygen_ranges), 0);
mean_pval = zeros(numel(oxygen_ranges), 0);
for i=1:numel(oxygen_ranges)-1
    mean_rho(i) = mean(diag(corr_matrices{i}), 'all');
    mean_pval(i) = mean(diag(p_vals{i}), 'all');
end

%% plot for comparing gradients at different oxy levels

newcolors = viridis(16);
fig = figure('Renderer', 'painters', 'Position', [10 10 400 300]);
% plot of relative error of J_ox measurement as a function of outside
% oxygen level

rel_err = zeros(numel(precise_oxygen_levels), 10);
for k=1:numel(oxygen_ranges)
    dist_all=oxygen_ranges_data{k}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{k}.jox_cell_kn_all;
    jox = mean(jox_cell_kn_all, 'omitnan');
    sigma_jox = std(jox_cell_kn_all, 'omitnan')./sqrt(size(jox_cell_kn_all,1));
    rel_err(k, :) = sigma_jox./jox;
end

subplot(2, 2, 1);
hold on;
plot(precise_oxygen_levels, max(rel_err, [], 2), 'o', 'LineWidth',1.5, ...
    'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('max(\sigma_{J_{ox}}/J_{ox})');
set(gca,'FontSize',12);
xlim([0 12])

% plot of correlation coefficient between data set at current
% with data set of next oxygen level
subplot(2, 2, 2);
hold on;
plot(precise_oxygen_levels(1:end-1), mean_rho, 'o', 'LineWidth',1.5,...
    'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('Kendalls tau correlation coefficient');
set(gca,'FontSize',12);
xlim([0 12])

% plot of mean squared distance to next higher J_ox data set
subplot(2, 2, 3);
hold on;
plot(precise_oxygen_levels(1:end-1), msds, 'o', 'LineWidth',1.5,...
    'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('mean squared distance to next set');
set(gca,'FontSize',12);
xlim([0 12])

% plot of KS test results (test, whether current and next J_ox data set
% where sampled from the same prob. distr.
subplot(2, 2, 4);
hold on;
plot(precise_oxygen_levels(1:end-1), mean(kstest2_results, 2),...
    'o', 'LineWidth',1.5, 'DisplayName','rejected h0','MarkerSize',10)

plot(precise_oxygen_levels(1:end-1), mean(kstest2_p, 2), 'o', 'LineWidth',1.5,...
    'Color','red', 'DisplayName','confidence level p', 'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('% of rejected h0 hypothesis (KS test)');
set(gca,'FontSize',12);
xlim([0 12])
legend()

savefig(string(plot_path)+'Jox_vs_dist_empComp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_empComp'+string(temp_string)+'.png')