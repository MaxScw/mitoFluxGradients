%% load setup file for file paths and data
clc
clear all
temp_resolved = false;
run("setup.m")

%% reshape raw data for easier handling

% choose number of rings used in spatial analysis
n_rings = 10;

if temp_resolved==true
diff_coeffs = diffusivity(temp_ind);%[1600, 1700, 1800, 1900, 2000, 2100, 2150, 2200, 2250, 2300, 2400, 2450, 2500, 2550, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500];
else
diff_coeffs = diffusivity(4);
end

% first, the data (J_ox, r) is reshaped into matrices of dimension
% (n_oxy_levels)x(n_rings) for easier access

jox_data = zeros(numel(precise_oxygen_levels), n_rings);
r_data = zeros(numel(precise_oxygen_levels), n_rings);
sigma_jox_data = zeros(numel(precise_oxygen_levels), n_rings);
sigma_r_data = zeros(numel(precise_oxygen_levels), n_rings);

jox_data_mitoDist = zeros(numel(precise_oxygen_levels), n_rings);
sigma_jox_data_mitoDist = zeros(numel(precise_oxygen_levels), n_rings);
mean_mito_density = mean(ratio_num_all_mitotracker, 1);%mean(density_num_all_mitotracker, 1);

stderr_mito_density = std(density_num_all_mitotracker, 1)./sqrt(numel(density_num_all_mitotracker(:, 1)));

for ind=1:numel(precise_oxygen_levels)
    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    jox_data(ind, :) = data(2, :);
    r_data(ind, :) = data(1, :);
    sigma_jox_data(ind, :) = data_stderr(2, :);
    sigma_r_data(ind, :) = data_stderr(1, :);

    % APPLY mito_weight TO data

    % include inheterogeneity
    mito_weight = mean_mito_density;
    
    % average over inhet.
    %mito_weight = mean(mean_mito_density);
    
    % apply conversion factor to obtain J_consum
    jox_data_mitoDist(ind, :) = data(2, :).*mito_weight./2;    

end 

%% load jox integration results

mitoCorr = true;
if mitoCorr==true
    corr_string = '_corrected';
else
    corr_string = '';
end

% load v_max, k_m profiles
load(pp_path + "v_max_profiles"+corr_string+string(temp_string)+".mat")
load(pp_path + "v_max_profiles"+corr_string+"_sigma"+string(temp_string)+".mat")
load(pp_path + "k_m_profiles"+corr_string+string(temp_string)+".mat")
load(pp_path + "k_m_profiles"+corr_string+"_sigma"+string(temp_string)+".mat")

% save corrected oxygen levels
load(pp_path + "cOxy"+corr_string+string(temp_string)+".mat")

% save predicted J_ox gradient
load(pp_path + "jox_pred"+string(temp_string)+".mat")


% save corrected oxygen levels with fitted spatial kinetics
load(pp_path + "coxy_pred_theory"+corr_string+string(temp_string)+".mat")

% save predicted J_ox gradient with fitted spatial kinetics
load(pp_path + "jox_pred_theory"+corr_string+string(temp_string)+".mat")

%% plot corrected oxygen levels


fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],'color','white');

cs = flip(viridis(20));
cs2 = magma(4);
colormap(cs);

% subplot(1, 2, 1)
hold on
for ind=[2, 3, 4, 5, 7, 10, 12, 14, 16]
    plot(r_data(ind, :), squeeze(cOxy_all(2, ind, :))./squeeze(cOxy_all(2, ind, end)), 'o',...
         'Color',cs(ind, :), 'MarkerSize',10, 'LineWidth',1.5)
end
xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex')
ylabel('$\hat{c}(c_\mathrm{out}, r)/\hat{c}(c_\mathrm{out}, R)$', 'Interpreter','latex')
set(gca,'FontSize',19)
title('norm. subcell. oxygen levels', 'Interpreter','latex')

cb = colorbar;
% clim([1 16])
% title(cb, '$c_\mathrm{out} (\mu M)$', 'Interpreter', 'latex')
% cb.Ticks = [2, 3, 4, 5, 7, 10, 12, 14, 16];
% cb.TickLabels = round(precise_oxygen_levels([1, 2, 3, 4, 5, 7, 10, 12, 14, 16]), 2);

clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$c_\mathrm{out}$ [\,\,\,\,M]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);

ylim([0.65, 1.01])
%yscale('log')
xlim([0, 36])

savefig(string(plot_path)+'coxy_vs_dist_int'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'coxy_vs_dist_int'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'coxy_vs_dist_int'+string(temp_string)+'.eps')


