%% load setup file for file paths and data
clc
clear all
temp_resolved = true;
run("setup.m")


%% plot highest oxy level jox gradient for all temps
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],...
    'Color','white');
cs = plasma(4+1);
cs = cs(1:end-1, 1:end);
colormap(cs)

hold on

xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex');
ylabel('$J_{\textrm{ox}}$ [$\mu$ M/s]', 'Interpreter','latex');

set(gca,'FontSize',19);

cb = colorbar;
clim([0.5 4.5])
title(cb, 'T [$^\circ$C]', 'Interpreter', 'latex')

cb.Ticks = linspace(1, 4, 4);
kelvin_temps = temperature(1:end-1)+273.15;
cb.TickLabels = temperature(1:end-1);

max_oxy = 50; % muM
min_oxy = 20;
title_string = 'average over '+string(min_oxy)+'$< c_\mathrm{out} <$'+string(max_oxy)+'$\mu$ M';
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
