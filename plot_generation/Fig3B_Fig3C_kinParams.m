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







%% plot v_max(r), k_m(r) profiles from data + approximate functional fits

fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
    'Color','white');

fs = 19;
cs = jet(16);

% v_max profile

hold on
% plot(data(1, 2:end), vMax_all(2, 2:end), 'HandleVisibility','off',...
%       'LineWidth',1.5, 'Color','red')
errorbar(data(1, 2:end), vMax_all(2, 2:end), sigma_vMax_all(2, 2:end),...
         sigma_vMax_all(2, 2:end),'o', 'MarkerSize',10,...
         'LineWidth',1.5, 'Color','red')

y_weights = 1./sigma_vMax_all(2, 1:end);

% set parameters to fit, increment and initial values (LAMBDA LINEAR MODEL)
% param = [vmax, km]

dp = [0.001 0.001 0.001];19

p0 = [65 0.0001 0];

pmin = [-200 1e-8 -200];
pmax = [200 200 200];

% fit

fit_start = 2;
fit_end = 10;
max_iter = 1000;

[p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                 lm(@exp_model, p0,...
                 data(1, 1:end)', vMax_all(2, 1:end)', y_weights(fit_start:fit_end), dp,...
                 pmin, pmax, [], fit_start, fit_end, max_iter);
r_range = linspace(data(1, 2), data(1, end), 100);
%legendstring = '$a*exp(r/b)+c$ (a='+string(round(p(1), 3))+',b='+string(round(1/p(2), 2))+',c='+string(round(p(3), 2))+')';
legendstring = '$a\cdot exp(r/b)+c$';
% plot(data(1, 2:end), exp_model(data(1, 2:end), p), 'LineWidth',1.5, ...
%      'Color','black', 'LineStyle','--')

p_vmax = p;

xlabel('distance from cell center ($r$) [\,\,\,\,m]','Interpreter','latex')
ylabel('$v_{\mathrm{max}}$ [\,\,\,\,M/s]', 'Interpreter','latex')

set(gca,'FontSize',fs)

% legend({'$v_{\mathrm{max}}$', legendstring}, 'Location','northwest', 'Interpreter','latex')
title('$v_{\mathrm{max}} (r)$ profile', 'Interpreter','latex')

savefig(string(plot_path)+'vMaxFit_vs_dist'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'vMaxFit_vs_dist'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'vMaxFit_vs_dist'+string(temp_string)+'.eps')


% k_m profile

fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
    'Color','white');


hold on
% plot(data(1, 2:end), kM_all(2, 2:end), 'HandleVisibility','off',...
%       'LineWidth',1.5, 'Color','red')
errorbar(data(1, 2:end), kM_all(2, 2:end), sigma_kM_all(2, 2:end),...
         sigma_kM_all(2, 2:end),'o', 'MarkerSize',10,...
         'DisplayName','$K_{\mathrm{M}}$',  'LineWidth',1.5, 'Color','red')

y_weights = 1./sigma_kM_all(2, 2:end);

% set parameters to fit, increment and initial values (LAMBDA LINEAR MODEL)
% param = [vmax, km]

dp = [0.001 0.001];

p0 = [0.001 0.1];

pmin = [-200 -200];
pmax = [200 200];

% fit

fit_start = 1;
fit_end = numel(sigma_vMax_all(1, 2:end));
max_iter = 1000;

[p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                 lm(@linear_model, p0,...
                 data(1, 2:end)', kM_all(2, 2:end)', y_weights, dp,...
                 pmin, pmax, [], fit_start, fit_end, max_iter);
r_range = linspace(data(1, 2), data(1, end), 100);
%legendstring = '$a*r+b$ (a='+string(round(p(1), 3))+',b='+string(round(p(2), 2))+')';
legendstring = '$a\cdot r+b$';
% plot(data(1, 2:end), linear_model(data(1, 2:end), p), 'LineWidth',1.5, ...
%      'Color','black','LineStyle', '--', 'DisplayName',legendstring)

xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex')
ylabel('$K_{\mathrm{M}}$ [\,\,\,\,M]', 'Interpreter','latex')

P_km = p;

set(gca,'FontSize',fs)
% legend('Location','best', 'Interpreter','latex')
title('$K_{\mathrm{M}}(r)$ profile', 'Interpreter','latex')

savefig(string(plot_path)+'KmFit_vs_dist'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'KmFit_vs_dist'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'KmFit_vs_dist'+string(temp_string)+'.eps')
