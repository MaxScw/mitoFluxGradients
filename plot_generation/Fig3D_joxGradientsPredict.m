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





%% plot J_ox(r) predicted from v_max(r), k_m(r)

% plot integrated jox for all oxygen levels
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],...
    'Color','white');
fs = 22;

start_ring = 2;
cs = flip(viridis(20));
colormap(cs)
X2_all = [];

%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
% plot all oxygen ranges together 
for oxy=1:16
    hold on
    % plot(r_data(oxy, :), squeeze(jox_pred(2, oxy, :)), 'Color',cs(oxy, :), ...
    %      'LineWidth',1.5)
    if oxy==10
        vis = 'on';
    else
        vis = 'off';
    end
    

    %joxPred = squeeze(cOxy_all(2, oxy, :)).*vMax_all(2, :)'./(squeeze(cOxy_all(2, oxy, :))+ kM_all(2, :)');
    %joxPred = squeeze(cOxy_all(2, oxy, :)).*exp_model(r_data(oxy, :), p_vmax)'./(squeeze(cOxy_all(2, oxy, :))+ linear_model(r_data(oxy, :)', P_km));
    plot(r_data(oxy, start_ring:end), squeeze(jox_int_allOxy(1, oxy, :)),...
         'Color',cs(oxy, :), 'LineWidth',1.5, 'HandleVisibility',vis, ...
         'DisplayName','$J_{\mathrm{ox}}^{\mathrm{theory}}$')

    errorbar(r_data(oxy, 2:end), jox_data(oxy, 2:end), sigma_jox_data(oxy, 2:end), ...
        sigma_jox_data(oxy, 2:end),'o', 'Color',cs(oxy, :), 'LineWidth',1.5, ...
        'MarkerSize', 10, 'HandleVisibility',vis, 'DisplayName','$J_{\mathrm{ox}}$')

    chi_sq = mean((jox_data(oxy, :) - squeeze(jox_int_allOxy(1, oxy, start_ring:end))).^2);

    X2_all = [X2_all chi_sq];
end

xlim([4, 36])
ylim([-25, 140])

xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex');
ylabel('ETC flux [$\mu$ M/s]', 'Interpreter','latex');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',fs);

cb = colorbar;

clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$c_\mathrm{out}$ [$\mu$ M]', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);
legend('Location','northwest', 'Interpreter','latex')

savefig(string(plot_path)+'Jox_vs_dist_int'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_int'+string(temp_string)+'.png')



% plot integrated jox for maximally different oxygen levels
fig = figure('Renderer', 'painters', 'Position', [10 10 600 400]);

start_ring = 2;
cs = flip(viridis(20));
colormap(cs)
X2 = [];

%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
% plot all oxygen ranges together 
for oxy=[1, 2, 16]
    hold on
    % plot(r_data(oxy, :), squeeze(jox_pred(2, oxy, :)), 'Color',cs(oxy, :), ...
    %      'LineWidth',1.5)
    errorbar(r_data(oxy, 2:end), jox_data(oxy, 2:end), sigma_jox_data(oxy, 2:end), ...
        sigma_jox_data(oxy, 2:end),'o', 'Color',cs(oxy, :), 'LineWidth',1.5, ...
        'MarkerSize', 10)

    %joxPred = squeeze(cOxy_all(2, oxy, :)).*vMax_all(2, :)'./(squeeze(cOxy_all(2, oxy, :))+ kM_all(2, :)');
    joxPred = squeeze(jox_int_allOxy(1, oxy, :));
    plot(r_data(oxy, start_ring:end), joxPred,...
         'Color',cs(oxy, :), 'LineWidth',1.5)

    chi_sq = sum(((jox_data(oxy, start_ring:end) - joxPred').^2)./jox_data(oxy, start_ring:end));

    X2 = [X2 chi_sq];
end

xlim([4, 36])
ylim([-25, 140])
mean(X2)
xlabel('distance to oocyte center ($\mu$ m)', 'Interpreter','latex');
ylabel('predicted $J_{ox}$ ($\mu$ M/s)', 'Interpreter','latex');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',fs);

cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, '$c_\mathrm{out} (\mu M)$', 'Interpreter', 'latex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);

savefig(string(plot_path)+'Jox_vs_dist_int_maxDiff'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_int_maxDiff'+string(temp_string)+'.png')
