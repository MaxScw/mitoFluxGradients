%% load setup file for file paths and data
clc
clear all
temp_resolved = true;
run("setup.m")

%% load in single-temp data
for temp_ind=1:4
    
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
    
    
    % plot NADH free and bound for different external oxygen for all rings
    
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
    
        xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex');
        ylabel('$[\mathrm{NADH}_\mathrm{f}]$ [$\mu$ M]', 'Interpreter','latex');
    
        
    
        set(gca,'FontSize',19);
    end
    
    if (temp_resolved==true)&&(temp_ind==1)
    cb = colorbar;
    clim([1.5 numel(precise_oxygen_levels)+0.5])
    title(cb, '$c_\mathrm{out}$ [$\mu$ M]', 'Interpreter', 'latex')
    cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
    cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);
    end
    
    if temp_resolved==true
    title('T='+string(temperature(temp_ind))+'$^\circ C$', 'Interpreter', 'latex')
    end
    
    xlim([4, 36])
    %ylim([20 85])
    
    
    savefig(string(plot_path)+'NADHf_vs_dist_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'NADHf_vs_dist_emp'+string(temp_string)+'.png')
    
    
    
    
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
    
        xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex');
        ylabel('$[\mathrm{NADH}_\mathrm{b}]$ [$\mu$ M]', 'Interpreter','latex');
    
        
    
        set(gca,'FontSize',19);
    end
    
    if (temp_resolved==true)&&(temp_ind==1)
    cb = colorbar;
    clim([1.5 numel(precise_oxygen_levels)+0.5])
    title(cb, '$\mathrm{c}_{\mathrm{out}}$ [$\mu$ M]', 'Interpreter', 'latex')
    cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
    cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);
    end
    
    if temp_resolved==true
    title('T='+string(temperature(temp_ind))+'$^\circ C$', 'Interpreter', 'latex')
    end
    
    xlim([4, 36])
    %ylim([15 40])
    
    savefig(string(plot_path)+'NADHb_vs_dist_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'NADHb_vs_dist_emp'+string(temp_string)+'.png')
    
end