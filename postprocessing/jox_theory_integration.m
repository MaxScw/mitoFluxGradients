%% designation of oxygen ranges to analyse
clear all;
clc;

oxygen_ranges={"run_oxy_l0_h0p1", "run_oxy_l0p1_h0p2", "run_oxy_l0p2_h0p3", "run_oxy_l0p3_h0p4",...
               "run_oxy_l0p4_h0p6", "run_oxy_l0p6_h0p8",...
               "run_oxy_l0p8_h1", "run_oxy_l1_h1p2",...
               "run_oxy_l1p2_h1p4", "run_oxy_l1p4_h1p6",...
               "run_oxy_l1p6_h1p8",...
               "run_oxy_l1p8_h2", "run_oxy_l2_h2p8",...
               "run_oxy_l2p8_h3p6", "run_oxy_l3p6_h4p4",...
               "run_oxy_l4p4_h5p2"};
num_oxygen_ranges={[0, 0.1], ...
                   [0.1, 0.2], [0.2, 0.3], [0.3, 0.4],...
                   [0.4, 0.6], [0.6, 0.8], [0.8, 1.0], [1.0, 1.2],...
                   [1.2, 1.4], [1.4, 1.6], [1.6, 1.8], [1.8, 2.0],...
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4],... 
                   [4.4, 5.2]};

temp_resolved = false;
data_path = "/home/mx/mitoFluxGradients/data/";
pp_path = "/home/mx/mitoFluxGradients/data/EXP_published/postprocessing_results/";
plot_path = "/home/mx/mitoFluxGradients/data/EXP_published/plots/";
% temp_resolved = true;
% data_path = "/home/mx/mitoFluxGradients/data/";
% pp_path = "../data/EXP_temperatureResolved/postprocessing_results/";
% plot_path = "../data/EXP_temperatureResolved/plots/";

%% loading of mito density, temperatures, solubilities, diffusivities

load(data_path + "/mito_density_oocytes/mito_density.mat")
mean_mito_density = mean(ratio_num_all_mitotracker, 1);

load(data_path + "exp_solubilities.mat")
load(data_path + "exp_temperatures.mat")
load(data_path + "exp_diffusivities.mat")

if temp_resolved==true
    temp_ind = 2;
    temp_string = "_T"+string(temperature(temp_ind))+"C"
else 
    temp_string = "";
end

%% loading of read-in data

load(pp_path + "plot_data_multiple_oxy_ranges"+temp_string+".mat");

precise_oxygen_levels = [];
sigma_precise_oxygen_levels = [];
for oxy_ind=1:numel(oxygen_ranges)
    oxygen_ranges_data{oxy_ind}.o2_levels;
    precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
    sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, mean(oxygen_ranges_data{oxy_ind}.sigma_o2_levels, "omitnan")];
end


if temp_resolved==true
precise_oxygen_levels = solubility(temp_ind).*precise_oxygen_levels./20.946;
sigma_precise_oxygen_levels = solubility(temp_ind).*sigma_precise_oxygen_levels./20.946;
sigma_precise_oxygen_levels(sigma_precise_oxygen_levels==0) = mean(sigma_precise_oxygen_levels);
else
precise_oxygen_levels = (213.5/20.946).*precise_oxygen_levels;
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

%% reshape raw data for easier handling

% choose number of rings used in spatial analysis
n_rings = 10;

if temp_resolved==true
diff_coeffs = diffusivity(temp_ind);%[1600, 1700, 1800, 1900, 2000, 2100, 2150, 2200, 2250, 2300, 2400, 2450, 2500, 2550, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500];
else
ref_diff_coeff = 3320;%2910;
diff_coeffs = ref_diff_coeff+[-400 -300 -200 -100 -50 0 50 100 200 300 400];
%diff_coeffs = [1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 2000 2250 2500 2750 3000 4000];
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

%% plot v_max(r), k_m(r) profiles from data + approximate functional fits
% 
% fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
%     'Color','white');
% cs = jet(16);

% v_max profile
% 
% hold on
% plot(data(1, 2:end), vMax_all(2, 2:end), 'HandleVisibility','off',...
%       'LineWidth',1.5, 'Color','red')
% errorbar(data(1, 2:end), vMax_all(2, 2:end), sigma_vMax_all(2, 2:end),...
%          sigma_vMax_all(2, 2:end),'o', 'MarkerSize',10,...
%          'LineWidth',1.5, 'Color','red')
% 
y_weights = 1./sigma_vMax_all(2, 1:end);

% set parameters to fit, increment and initial values (LAMBDA LINEAR MODEL)
%param = [vmax, km]

dp = [0.001 0.001 0.001];

p0 = [65 0.0001 0];

pmin = [-200 1e-8 -200];
pmax = [200 200 200];

% fit

fit_start = 1;
fit_end = 10;
max_iter = 1000;

[p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                 lm(@exp_model, p0,...
                 data(1, 1:end)', vMax_all(2, 1:end)', y_weights(fit_start:fit_end), dp,...
                 pmin, pmax, [], fit_start, fit_end, max_iter);
% r_range = linspace(data(1, 2), data(1, end), 100);
% %legendstring = '$a*exp(r/b)+c$ (a='+string(round(p(1), 3))+',b='+string(round(1/p(2), 2))+',c='+string(round(p(3), 2))+')';
% legendstring = '$a\cdot exp(r/b)+c$, $R^2=$'+ string(round(R_sq, 2));
% plot(data(1, 2:end), exp_model(data(1, 2:end), p), 'LineWidth',1.5, ...
%      'Color','black', 'LineStyle','--')
% 
p_vmax = p;
% 
% xlabel('distance from cell center ($r$) $[\mu m]$','Interpreter','latex')
% ylabel('$v_{\mathrm{max}} [\mu M/s]$', 'Interpreter','latex')
% 
% set(gca,'FontSize',19)
% 
% legend({'$v_{\mathrm{max}}$', legendstring}, 'Location','northwest', 'Interpreter','latex')
% title('$v_{\mathrm{max}}$ profile from mm fit', 'Interpreter','latex')


% k_m profile

% fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
%     'Color','white');
% 
% 
% hold on
% plot(data(1, 2:end), kM_all(2, 2:end), 'HandleVisibility','off',...
%       'LineWidth',1.5, 'Color','red')
% errorbar(data(1, 2:end), kM_all(2, 2:end), sigma_kM_all(2, 2:end),...
%          sigma_kM_all(2, 2:end),'o', 'MarkerSize',10,...
%          'DisplayName','$K_{\mathrm{M}}$',  'LineWidth',1.5, 'Color','red')
% 
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
% r_range = linspace(data(1, 2), data(1, end), 100);
% %legendstring = '$a*r+b$ (a='+string(round(p(1), 3))+',b='+string(round(p(2), 2))+')';
% legendstring = '$a\cdot r+b$, $R^2=$'+ string(round(R_sq, 2));
% plot(data(1, 2:end), linear_model(data(1, 2:end), p), 'LineWidth',1.5, ...
%      'Color','black','LineStyle', '--', 'DisplayName',legendstring)
% 
% xlabel('distance from cell center ($r$) $[\mu m]$', 'Interpreter','latex')
% ylabel('$K_{\mathrm{M}} [\mu M]$', 'Interpreter','latex')
% 
P_km = p;
% 
% set(gca,'FontSize',19)
% legend('Location','best', 'Interpreter','latex')
% title('$K_{\mathrm{M}}$ profile from mm fit', 'Interpreter','latex')

%% integrate RD equation with fitted params for spatial kinetics
first_ring = 1;

coxy_int_allOxy = zeros(numel(diff_coeffs), 16, n_rings-first_ring+1);
jox_int_allOxy = zeros(numel(diff_coeffs), 16, n_rings-first_ring+1);

all_x2 = zeros(16, numel(diff_coeffs));

corr_fact = (mean_mito_density(first_ring:end))./2;
param_profiles = false;
for diff_ind=1:numel(diff_coeffs)
    for ind=1:16
        D = diff_coeffs(diff_ind);
        c_star = precise_oxygen_levels(ind);
        params = [c_star, D];
        
        
        
        
        if param_profiles==true
            const = [corr_fact.*vMax_all(2, first_ring:end); kM_all(2, first_ring:end)];
        else
            const = [p_vmax; P_km; corr_fact'];
            
        end
        int_allOxy = theory_flux_integrator(data(1, first_ring:end), params, const, param_profiles);
        coxy_int_allOxy(diff_ind, ind, :) = int_allOxy(1, :);

        if param_profiles==false
            jox_int_allOxy(diff_ind, ind, :) = int_allOxy(1, :).*exp_model(data(1, first_ring:end), p_vmax)./(int_allOxy(1, :) + linear_model(data(1, first_ring:end), P_km));
        else
            jox_int_allOxy(diff_ind, ind, :) = int_allOxy(1, :).*vMax_all(2, first_ring:end)./(int_allOxy(1, :) + kM_all(2, first_ring:end));
        end
    
        all_x2(ind, diff_ind) = sum(((jox_data(ind, first_ring:end)...
                                - squeeze(jox_int_allOxy(diff_ind, ind, :))').^2)...
                                ./jox_data(ind, first_ring:end));
    end
end
%%
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],'color','white');
set(gca,'FontSize',19)
hold on
yyaxis left
plot(diff_coeffs, mean(all_x2, 1), 'LineWidth',1.5, 'LineStyle',':', 'DisplayName','spatial RD sol.')
xline(diffusivity(4), 'LineStyle','--', 'Color','black', 'LineWidth',1.5, 'DisplayName','$D=3320 \mu m^2/s$')
ylabel('$\overline{\chi^2}$', 'Interpreter','latex')

yyaxis right
coxy_var = mean((coxy_int_allOxy(2:end, :, :) - coxy_int_allOxy(1:end-1, :, :)), [2, 3])./(diff_coeffs(2:end)-diff_coeffs(1:end-1))';
plot(diff_coeffs(2:end), coxy_var, 'Color','red', 'LineWidth',1.5, 'HandleVisibility','off')
ylabel('$\Delta \hat{c}/\Delta D$ [$\mu M/(\mu m^2/s)$]', 'Interpreter','latex')

xlabel('D [$\mu m^2/s$]', 'Interpreter','latex')

xlim([100, diff_coeffs(end)])
legend('Interpreter','latex', Location='east')

savefig(string(plot_path)+'chiSqSol_Drobust'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'chiSqSol_Drobust'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'chiSqSol_Drobust'+string(temp_string)+'.eps')
%%
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],'color','white');
hold on
cs = summer(numel(diff_coeffs));
oxy_ind = 2;
precise_oxygen_levels(oxy_ind)
count = 0;
for ind=1:numel(diff_coeffs)
    count = count + 1;
plot(r_data(1, first_ring:end), squeeze(coxy_int_allOxy(ind, oxy_ind, :))./squeeze(coxy_int_allOxy(ind, oxy_ind, end)),'o', 'color', cs(count, :),...
    'MarkerSize',10, 'LineWidth',1.5)
end

xlabel('distance from cell center ($r$) [\,\,\,\,m]', 'Interpreter','latex')
ylabel('$\hat{c}(c_\mathrm{out}, r)/\hat{c}(c_\mathrm{out}, R)$', 'Interpreter','latex')
set(gca,'FontSize',19)

colormap(cs);
cb = colorbar;
clim([diff_coeffs(1) diff_coeffs(end)])
title(cb, '$D$ [\,\,\,\,m$^2$/s]', 'Interpreter', 'latex')
%cb.Ticks = [diff_coeffs(1) diff_coeffs(5) diff_coeffs(10) diff_coeffs(15) diff_coeffs(22)];
%cb.TickLabels = [diff_coeffs(1) diff_coeffs(5) diff_coeffs(10) diff_coeffs(15) diff_coeffs(22)];

ylim([0.65, 1.01])
xlim([0, 36])

level = ref_diff_coeff;
h_axes = axes('position', cb.Position, 'ylim', cb.Limits, 'color', 'none', 'visible','off');
lim = [1 10];
line(lim, level*[1 1], 'color', 'red', 'parent', h_axes, 'Linewidth', 5);

savefig(string(plot_path)+'coxy_vs_dist_int_Drobust'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'coxy_vs_dist_int_Drobust'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'coxy_vs_dist_int_Drobust'+string(temp_string)+'.eps')

%% save resulting J_ox^theory profiles

% save predicted J_ox gradient
save(pp_path + "coxy_pred_theory"+corr_string+string(temp_string)+".mat", "coxy_int_allOxy")
save(pp_path + "jox_pred_theory"+corr_string+string(temp_string)+".mat", "jox_int_allOxy")

%% obtain decay lengths from oxygen gradient by fitting sinh(R/lambda)
figure(2)
hold on
cs = viridis(16);

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];
diff_ind = 1;

for oxy=1:numel(precise_oxygen_levels)

    % p = [R, A, lambda]
    p0 = [r_data(oxy, end) 1 30];
    dp = [0 0.001 0.001];
    
    y_weights = squeeze(coxy_int_allOxy(diff_ind, oxy, :))./sigma_precise_oxygen_levels;
    fit_start = 2;
    fit_end = 10;

    pmin = [0 1e-4 1e-4];
    pmax = [1e4 1e4 1e4];

    c = [];

    % fit
    max_iter = 100;
   
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@rd_solver_linear, p0,...
                     r_data(oxy, :)', squeeze(coxy_int_allOxy(diff_ind, oxy, :))./squeeze(coxy_int_allOxy(diff_ind, oxy, end)), y_weights(fit_start:fit_end), dp,...
                     pmin, pmax, c, fit_start, fit_end, max_iter);
    % !!!optional:
    % use sum of squared residuals instead of built-in method of lm to evaluate X2
    %x2_jox = sum(((r_data(oxy,fit_start:fit_end) - linear_model(norm_log_coxy(fit_start:fit_end), p, c)).^2)./norm_log_coxy(fit_start:fit_end));
    %X2 = x2_jox;
    
    if oxy<4
    plot(r_data(oxy, :), squeeze(coxy_int_allOxy(diff_ind, oxy, :))./squeeze(coxy_int_allOxy(diff_ind, oxy, end)), "o", Color=cs(oxy, :))
    plot(r_data(oxy, :), rd_solver_linear(r_data(oxy, :), p, c), Color=cs(oxy, :), LineWidth=1.5)
    end

    fit_params = [fit_params, p];
    sigma_fit_params = [sigma_fit_params, sigma_p];
    chi_squareds = [chi_squareds X2];
    R_squareds = [R_squareds R_sq];
end

decay_lengths = fit_params(3, :)
sigma_decay_lengths = decay_lengths.*(sigma_fit_params(3, :)./fit_params(3, :));

%save(pp_path + "sinhDecayLength_coxy" + temp_string + ".mat", "decay_lengths")
%save(pp_path + "sigma_sinhDecayLength_coxy" + temp_string + ".mat", "sigma_decay_lengths")
