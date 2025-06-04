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
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4],... 
                   [4.4, 5.2]};

temp_resolved = false;
data_path = "/home/mx/mitoFluxGradients/data/";
pp_path = '../data/EXP_published/postprocessing_results/';
plot_path = '../data/EXP_published/plots/';
% temp_resolved = true;
% data_path = "/home/mx/mitoFluxGradients/data/";
% pp_path = '../data/EXP_temperatureResolved/postprocessing_results/';
% plot_path = '../data/EXP_temperatureResolved/plots/';

%% loading of mito density, temperatures, solubilities, diffusivities
load(data_path + "mito_density_oocytes/mito_density.mat")
mean_mito_density = mean(ratio_num_all_mitotracker, 1);

load(data_path + "exp_solubilities.mat")
load(data_path + "exp_temperatures.mat")
load(data_path + "exp_diffusivities.mat")

if temp_resolved==true
    temp_ind = 4;
    temp_string = '_T'+string(temperature(temp_ind))+'C';
else 
    temp_string = '';
end
%% loading of read-in (experimental) data

% choose which temperature to plot analysis for (from T=[22, 28, 31, 36]C)
load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');

precise_oxygen_levels = [];
sigma_precise_oxygen_levels = [];
for oxy_ind=1:numel(oxygen_ranges)
    precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
    sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, mean(oxygen_ranges_data{oxy_ind}.sigma_o2_levels, 'omitnan')];
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

%% load results from J_ox integration

% load v_max, k_m profiles
load(string(pp_path)+'v_max_profiles_corrected'+string(temp_string)+'.mat')
load(string(pp_path)+'v_max_profiles_corrected_sigma'+string(temp_string)+'.mat')
load(string(pp_path)+'k_m_profiles_corrected'+string(temp_string)+'.mat')
load(string(pp_path)+'k_m_profiles_corrected_sigma'+string(temp_string)+'.mat')

% load corrected oxygen levels
load(string(pp_path)+'cOxy_corr'+string(temp_string)+'.mat')

% load predicted J_ox gradient
load(string(pp_path)+'jox_pred'+string(temp_string)+'.mat')




%% plot corrected oxygen levels

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
cs = viridis(16);
cs2 = magma(4);
colormap(cs);

subplot(1, 2, 1)
hold on
for ind=[1, 2, 3, 4, 5, 7, 10, 12, 14, 16]
    plot(r_data(ind, :), squeeze(cOxy_all(2, ind, :))./squeeze(cOxy_all(2, ind, end)), 'o',...
         'Color',cs(ind, :), 'MarkerSize',10, 'LineWidth',1.5)
end
xlabel('distance from cell center (\mu m)')
ylabel('c(c^*, r)/c(c^*, R)')
set(gca,'FontSize',15)
title('norm. subcell. oxygen levels')

cb = colorbar;
clim([1 16])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')
cb.Ticks = [1, 2, 3, 4, 5, 7, 10, 12, 14, 16];
cb.TickLabels = round(precise_oxygen_levels([1, 2, 3, 4, 5, 7, 10, 12, 14, 16]), 2);
ylim([0.65, 1.01])
%yscale('log')
xlim([0, 36])

subplot(1, 2, 2)
plot(precise_oxygen_levels, mean(abs((squeeze(cOxy_all(2, :, :))-precise_oxygen_levels')), 2)./precise_oxygen_levels',...
     'o', 'MarkerSize',10, 'LineWidth',1.5)
ylabel('\langle |c(c^*, r) - c^*| \rangle_r /c^*')
xlabel('external oxygen conc. c* (\mu M)')
set(gca,'FontSize',15)
title('rel. deviation from flat profile c(r)=c^*')

savefig(string(plot_path)+'coxy_vs_dist_int'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'coxy_vs_dist_int'+string(temp_string)+'.png')

%% Analyse J_ox/c_oxy as a function of cell center distance

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
ax = gca; 
ax.FontSize = 20;
cs = viridis(16);
colormap(cs)
hold on

start_oxy = 2;
end_oxy =16;

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];

for oxy=start_oxy:end_oxy

    jox_over_coxy = jox_data(oxy, :)./squeeze(cOxy_all(2, oxy, :))';
    sigma_jox_over_coxy = sigma_jox_data(oxy, :)./squeeze(cOxy_all(2, oxy, :))';
    
    if oxy < 17
        subplot(1, 2, 1)
        hold on
        title('\gamma_{eff}(r, c^*)')
    else
        subplot(1, 2, 2)
        hold on
        title('high outside oxygen')
    end
  
    plot(r_data(oxy, :), jox_over_coxy,...
         'Color',cs(oxy, :), 'LineWidth',1.5, 'MarkerSize',15)
    xlabel('dist. from oocyte center (\mu m)')
    ylabel('J_{ox}(r) / c(r) (1/s)')
    %title('J_{ox}(r)/c(r)')%, \rho_{mito} corr.')
    %ylim([-50, 250])
    p0 = [0 10];
    dp = [0 0.001];
    
    y_weights = 1./sigma_jox_over_coxy;
    fit_start = 2;
    fit_end = 10;
    
    pmin = [-1 -100];
    pmax = [1 100];

    c = [];

    % fit
    max_iter = 100;
    
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@linear_model, p0,...
                     r_data(oxy, :)', jox_over_coxy', y_weights(fit_start:fit_end), dp,...
                     pmin, pmax, c, fit_start, fit_end, max_iter);
    % !!!optional:
    % use sum of squared residuals instead of built-in method of lm to evaluate X2
    x2_jox = sum(((jox_over_coxy(fit_start:fit_end) - linear_model(r_data(oxy,fit_start:fit_end), p, c)).^2)./jox_over_coxy(fit_start:fit_end));
    X2 = x2_jox;

    fit_params = [fit_params, p];
    sigma_fit_params = [sigma_fit_params, sigma_p];
    chi_squareds = [chi_squareds X2];
    R_squareds = [R_squareds R_sq];

    hold on
    plot(r_data(oxy, :), linear_model(r_data(oxy, :), p), ...
         'Color',cs(oxy, :), 'LineStyle',':', 'LineWidth',1.5)
    % if oxy<3
    %     text(r_data(oxy, end),linear_model(r_data(oxy, end), p),string(oxy),'Color','red');
    %     text(r_data(oxy, 1)-1,jox_over_coxy(1),string(oxy),'Color','red');
    % end
end

cb = colorbar;
clim([1 16])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')
cb.Ticks = linspace(1, 16, 16);
cb.TickLabels = round(precise_oxygen_levels, 2);

%fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
subplot(1, 2, 2)
yyaxis left
set(gca,'YColor', 'blue');
plot(precise_oxygen_levels(start_oxy:end_oxy), fit_params(2, :), 'o', ...
     'MarkerSize',15, 'LineWidth',1.5, 'Color',	'blue', 'DisplayName','fit to J_{ox}(r)/c_{oxy}(r)')
hold on 
% plot(precise_oxygen_levels(start_oxy:end_oxy), vMax_all(2, end)./(kM_all(2, end)+precise_oxygen_levels(start_oxy:end_oxy)), 'o', ...
%      'MarkerSize',15, 'LineWidth',1.5, 'Color',	'green', 'DisplayName','v_{max}\cdot c_{oxy}(r=R)/k_m')
% legend()
ylabel('offset (1/s)')

yyaxis right
set(gca,'YColor','red');
bar(precise_oxygen_levels(start_oxy:end_oxy), chi_squareds, ...
    'FaceAlpha',0.3, 'EdgeAlpha',0, 'LineWidth',1.5, 'FaceColor','red', ...
    'HandleVisibility','off')
ylabel('\chi^2')
xlabel('outside oxygen level (\mu M)')
title('fit offset to J_{ox}(r)/c(r)')%, \rho_{mito} corr.')

savefig(string(plot_path)+'JoxOverCoxy_vs_dist_int'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'JoxOverCoxy_vs_dist_int'+string(temp_string)+'.png')

%% plot v_max(r), k_m(r) profiles from data + approximate functional fits

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
cs = jet(16);

% v_max profile
subplot(1, 2, 1)
hold on
plot(data(1, 2:end), vMax_all(2, 2:end), 'HandleVisibility','off',...
      'LineWidth',1.5, 'Color','blue')
errorbar(data(1, 2:end), vMax_all(2, 2:end), sigma_vMax_all(2, 2:end),...
         sigma_vMax_all(2, 2:end),'o', 'MarkerSize',10,...
         'LineWidth',1.5, 'Color','blue')

y_weights = 1./sigma_vMax_all(2, 1:end);

% set parameters to fit, increment and initial values (LAMBDA LINEAR MODEL)
% param = [vmax, km]

dp = [0.001 0.001 0.001];

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
legendstring = 'a*exp(r/b)+c (a='+string(round(p(1), 3))+',b='+string(round(1/p(2), 2))+',c='+string(round(p(3), 2))+')';
plot(data(1, 2:end), exp_model(data(1, 2:end), p), 'LineWidth',1.5, ...
     'Color','red')

p_vmax = p;

xlabel('distance from cell center (\mu m)')
ylabel('v_{max} (\mu M/s)')

set(gca,'FontSize',15)

legend({'v_{max}', legendstring}, 'Location','northwest')
title('v_{max} profile from mm fit')

% k_m profile
subplot(1, 2, 2)
hold on
plot(data(1, 2:end), kM_all(2, 2:end), 'HandleVisibility','off',...
      'LineWidth',1.5, 'Color','blue')
errorbar(data(1, 2:end), kM_all(2, 2:end), sigma_kM_all(2, 2:end),...
         sigma_kM_all(2, 2:end),'o', 'MarkerSize',10,...
         'DisplayName','k_{m}',  'LineWidth',1.5, 'Color','blue')

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
legendstring = 'a*r+b (a='+string(round(p(1), 3))+',b='+string(round(p(2), 2))+')';
plot(data(1, 2:end), linear_model(data(1, 2:end), p), 'LineWidth',1.5, ...
     'Color','red', 'DisplayName',legendstring)

xlabel('distance from cell center (\mu m)')
ylabel('k_{m} (\mu M)')

P_km = p;

set(gca,'FontSize',15)
legend('Location','best')
title('k_{m} profile from mm fit')

savefig(string(plot_path)+'vMaxKmFit_vs_dist'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'vMaxKmFit_vs_dist'+string(temp_string)+'.png')

%% plot correction in v_max(r), k_m(r) profiles

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
cs = jet(16);

subplot(1, 2, 1)
hold on
plot(data(1, 2:end), vMax_all(1, 2:end), 'Color','red','HandleVisibility','off',...
      'LineWidth',1.5)
errorbar(data(1, 2:end), vMax_all(1, 2:end), sigma_vMax_all(1, 2:end),...
         sigma_vMax_all(1, 2:end),'o','Color','red', 'MarkerSize',10,...
         'DisplayName','v_{max, 0}',  'LineWidth',1.5)
xlabel('distance from oocyte center (\mu m)')
ylabel('k_m (\mu M)')

plot(data(1, 2:end), vMax_all(2, 2:end), 'Color','green','HandleVisibility','off',...
     'LineWidth',1.5)
errorbar(data(1, 2:end), vMax_all(2, 2:end), sigma_vMax_all(2, 2:end),...
         sigma_vMax_all(2, 2:end),'o','Color','green', 'MarkerSize',10,...
         'DisplayName','v_{max, corr}',  'LineWidth',1.5)

set(gca,'FontSize',15)
legend('Location','best')
title('correction of v_{max} profiles by J_{ox} int.')

subplot(1, 2, 2)
hold on
plot(data(1, 2:end), kM_all(1, 2:end), 'Color','red','HandleVisibility','off',...
      'LineWidth',1.5)
errorbar(data(1, 2:end), kM_all(1, 2:end), sigma_kM_all(1, 2:end),...
         sigma_kM_all(1, 2:end),'o','Color','red', 'MarkerSize',10,...
         'DisplayName','k_{m, 0}',  'LineWidth',1.5)
xlabel('distance from oocyte center (\mu m)')
ylabel('v_{max} (1/s)')

plot(data(1, 2:end), kM_all(2, 2:end), 'Color','green','HandleVisibility','off',...
     'LineWidth',1.5)
errorbar(data(1, 2:end), kM_all(2, 2:end), sigma_kM_all(2, 2:end),...
         sigma_kM_all(2, 2:end),'o','Color','green', 'MarkerSize',10,...
         'DisplayName','k_{m, corr}', 'LineWidth',1.5)

set(gca,'FontSize',15)
legend('Location','best')
title('correction of k_{m} profiles by J_{ox} int.')

savefig(string(plot_path)+'vMaxKmCorr_vs_dist'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'vMaxKmCorr_vs_dist'+string(temp_string)+'.png')

%% plot J_ox(r) predicted from v_max(r), k_m(r)

% plot integrated jox for all oxygen levels
fig = figure('Renderer', 'painters', 'Position', [10 10 600 400]);

start_ring = 2;
cs = flip(viridis(20));
colormap(cs)
X2 = [];

%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
% plot all oxygen ranges together 
for oxy=1:16
    hold on
    % plot(r_data(oxy, :), squeeze(jox_pred(2, oxy, :)), 'Color',cs(oxy, :), ...
    %      'LineWidth',1.5)
    errorbar(r_data(oxy, 2:end), jox_data(oxy, 2:end), sigma_jox_data(oxy, 2:end), ...
        sigma_jox_data(oxy, 2:end),'o', 'Color',cs(oxy, :), 'LineWidth',1.5, ...
        'MarkerSize', 10)

    %joxPred = squeeze(cOxy_all(2, oxy, :)).*vMax_all(2, :)'./(squeeze(cOxy_all(2, oxy, :))+ kM_all(2, :)');
    joxPred = squeeze(cOxy_all(2, oxy, :)).*exp_model(r_data(oxy, :), p_vmax)'./(squeeze(cOxy_all(2, oxy, :))+ linear_model(r_data(oxy, :)', P_km));
    plot(r_data(oxy, start_ring:end), joxPred(start_ring:end),...
         'Color',cs(oxy, :), 'LineWidth',1.5)

    chi_sq = sum(((jox_data(oxy, :) - joxPred').^2)./jox_data(oxy, :));

    X2 = [X2 chi_sq];
end

xlim([4, 36])
ylim([-25, 140])

xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',15);

cb = colorbar;
clim([1 16])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 16, 16);
cb.TickLabels = round(precise_oxygen_levels, 2);

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
    joxPred = squeeze(cOxy_all(2, oxy, :)).*exp_model(r_data(oxy, :), p_vmax)'./(squeeze(cOxy_all(2, oxy, :))+ linear_model(r_data(oxy, :)', P_km));
    plot(r_data(oxy, start_ring:end), joxPred(start_ring:end),...
         'Color',cs(oxy, :), 'LineWidth',1.5)

    chi_sq = sum(((jox_data(oxy, :) - joxPred').^2)./jox_data(oxy, :));

    X2 = [X2 chi_sq];
end

xlim([4, 36])
ylim([-25, 140])
mean(X2)
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
%title('predicted J_{ox}(v_{max}(r), k_{m}(r), c(r))')
set(gca,'FontSize',15);

cb = colorbar;
clim([1 16])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 16, 16);
cb.TickLabels = round(precise_oxygen_levels, 2);

savefig(string(plot_path)+'Jox_vs_dist_int_maxDiff'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_int_maxDiff'+string(temp_string)+'.png')
%% plot predicted J_ox separated in different oxgen regimes (low, mid, high)

fig = figure('Renderer', 'painters', 'Position', [10 10 1800 400]);
cs = viridis(16);
colormap(cs)

%low
X2 = [];
for oxy=1:6
    subplot(1, 3, 1)
    hold on
    errorbar(r_data(oxy, 2:end), jox_data(oxy, 2:end), sigma_jox_data(oxy, 2:end), ...
        sigma_jox_data(oxy, 2:end),'o', 'Color',cs(oxy, :), 'LineWidth',1.5, ...
        'MarkerSize', 10)
    joxPred = squeeze(cOxy_all(2, oxy, :)).*exp_model(r_data(oxy, :), p_vmax)'./(squeeze(cOxy_all(2, oxy, :))+ linear_model(r_data(oxy, :)', P_km));
    plot(r_data(oxy, start_ring:end), joxPred(start_ring:end),...
         'Color',cs(oxy, :), 'LineWidth',1.5)

    chi_sq = sum(((jox_data(oxy, :) - joxPred').^2)./jox_data(oxy, :));
    X2 = [X2 chi_sq];
end
mean(X2)
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
set(gca,'FontSize',15);

%mid
X2 = [];
for oxy=7:12
    subplot(1, 3, 2)
    hold on
    errorbar(r_data(oxy, 2:end), jox_data(oxy, 2:end), sigma_jox_data(oxy, 2:end), ...
        sigma_jox_data(oxy, 2:end),'o', 'Color',cs(oxy, :), 'LineWidth',1.5, ...
        'MarkerSize', 10)
    joxPred = squeeze(cOxy_all(2, oxy, :)).*exp_model(r_data(oxy, :), p_vmax)'./(squeeze(cOxy_all(2, oxy, :))+ linear_model(r_data(oxy, :)', P_km));
    plot(r_data(oxy, start_ring:end), joxPred(start_ring:end),...
         'Color',cs(oxy, :), 'LineWidth',1.5)

    chi_sq = sum(((jox_data(oxy, :) - joxPred').^2)./jox_data(oxy, :));
    X2 = [X2 chi_sq];
end
mean(X2)
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
set(gca,'FontSize',15);

%high
X2 = [];
for oxy=13:16
    subplot(1, 3, 3)
    hold on
    errorbar(r_data(oxy, 2:end), jox_data(oxy, 2:end), sigma_jox_data(oxy, 2:end), ...
        sigma_jox_data(oxy, 2:end),'o', 'Color',cs(oxy, :), 'LineWidth',1.5, ...
        'MarkerSize', 10)
    joxPred = squeeze(cOxy_all(2, oxy, :)).*exp_model(r_data(oxy, :), p_vmax)'./(squeeze(cOxy_all(2, oxy, :))+ linear_model(r_data(oxy, :)', P_km));
    plot(r_data(oxy, start_ring:end), joxPred(start_ring:end),...
         'Color',cs(oxy, :), 'LineWidth',1.5)

    chi_sq = sum(((jox_data(oxy, :) - joxPred').^2)./jox_data(oxy, :));
    X2 = [X2 chi_sq];
end

mean(X2)
xlabel('distance to oocyte center (\mu m)');
ylabel('predicted J_{ox} (\mu M/s)');
set(gca,'FontSize',15);

cb = colorbar;
clim([precise_oxygen_levels(1) precise_oxygen_levels(end)])
title(cb, 'c_{oxy} (\mu M)', 'Interpreter', 'tex')
cb.Ticks = round(linspace(precise_oxygen_levels(1), precise_oxygen_levels(end), 10), 3);

savefig(string(plot_path)+'Jox_vs_dist_sep'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_sep'+string(temp_string)+'.png')

%% load all decay lengths obtained from linear RD (J_ox sinh fit), 
%  J_ox log-lin fit, c_oxy sinh fit, c_oxy log-lin fit

% load J_ox decay length obtained by sinh fit
linRD_fit = load(string(pp_path)+'gradient_fit_linear_model'+string(temp_string)+'.mat');
linRD_fit = linRD_fit.("fit_data_gradient");
linRD_fit_params = linRD_fit{1};
linRD_fit_sigma_params = linRD_fit{2};

linRD_dec = linRD_fit_params(3, :);
sigma_linRD_dec = linRD_fit_sigma_params(3, :);

% load J_ox decay length obtained by sinh fit
jox_fit = load(string(pp_path)+'sinhDecayLength_jox'+string(temp_string)+'.mat');
jox_sigmafit = load(string(pp_path)+'sigma_sinhDecayLength_jox'+string(temp_string)+'.mat');
jox_sinhDec = jox_fit.("decay_lengths");
sigma_jox_sinhDec = jox_sigmafit.("sigma_decay_lengths");

% load J_ox decay length obtained by log-lin fit
joxEmp_fit = load(string(pp_path)+'decayLength_joxEmp'+string(temp_string)+'.mat');
joxEmp_sigmafit = load(string(pp_path)+'sigma_decayLength_joxEmp'+string(temp_string)+'.mat');
joxEmp_dec = joxEmp_fit.("decay_lengths");
sigma_joxEmp_dec = joxEmp_sigmafit.("sigma_decay_lengths");

% load c_oxy decay length obtained by log-lin fit
% coxy_fit = load("decayLength_coxy.mat");
% coxy_sigmafit = load('sigma_decayLength_coxy.mat');
% coxy_dec = coxy_fit.("decay_lengths");
% sigma_coxy_dec = coxy_sigmafit.("sigma_decay_lengths");

% load c_oxy decay length obtained by sinh fit
coxy_sinhFit = load(string(pp_path)+'sinhDecayLength_coxy'+string(temp_string)+'.mat');
coxy_sigmaSinhFit = load(string(pp_path)+'sigma_sinhDecayLength_coxy'+string(temp_string)+'.mat');
coxy_sinhDec = coxy_sinhFit.("decay_lengths");
sigma_coxy_sinhDec = coxy_sigmaSinhFit.("sigma_decay_lengths");


%% plot all decay lengths (J_ox and c_oxy, both methods) as a function of oxygen
fig = figure(5);
hold on

% plot J_ox decay length obtained by sinh fit
errorbar(precise_oxygen_levels, jox_sinhDec, sigma_jox_sinhDec, sigma_jox_sinhDec,...
         'o', 'Color','green', ...
         'MarkerSize',15, 'LineWidth',1.5, 'DisplayName','J_{ox} sinh fit')

% plot J_ox decay length obtained by log_lin fit
errorbar(precise_oxygen_levels, joxEmp_dec, sigma_joxEmp_dec, sigma_joxEmp_dec,...
         'o', 'Color','red', ...
         'MarkerSize',15, 'LineWidth',1.5, 'DisplayName','J_{ox} direct fit')

% plot c_oxy decay length obtained by sinh fit
errorbar(precise_oxygen_levels, coxy_sinhDec, sigma_coxy_sinhDec, sigma_coxy_sinhDec,...
         'o', 'Color','cyan', ...
         'MarkerSize',15, 'LineWidth',1.5, 'DisplayName','c_{oxy} sinh fit')

% plot c_oxy decay length obtained by log-lin fit
% errorbar(precise_oxygen_levels, coxy_dec, sigma_coxy_dec, sigma_coxy_dec,...
%          'o', 'Color','blue', ...
%          'MarkerSize',15, 'LineWidth',1.5, 'DisplayName','c_{oxy} direct fit')

% plot some reference lines
plot(precise_oxygen_levels, 6*precise_oxygen_levels.^0.5, ...
     'LineWidth',1.5, 'LineStyle','--', 'Color','black', 'DisplayName','~c_{oxy, ext}^{0.5}')
plot(precise_oxygen_levels, 6*precise_oxygen_levels.^1, ...
     'LineWidth',1.5, 'LineStyle','--', 'Color','black', 'DisplayName','~c_{oxy, ext}^{1}')

legend('Location','southeast')
xlabel('c_{oxy, ext.} (\mu M)')
ylabel('\lambda(c_{oxy, ext})')
title('comparison of decay lengths \lambda', 'Interpreter','tex')

savefig(string(plot_path)+'Lambda_vs_coxy'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Lambda_vs_coxy'+string(temp_string)+'.png')
%%
plot(jox_sinhDec)

%% plot lambda(J_ox) as a function of lambda(c_oxy) for sinh fit method
fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
set(gca,'FontSize',15)
hold on

cs = flip(viridis(20));
colormap(cs)

for oxy=1:16
errorbar(coxy_sinhDec(oxy), jox_sinhDec(oxy), sigma_coxy_sinhDec(oxy), sigma_coxy_sinhDec(oxy),...
    sigma_jox_sinhDec(oxy), sigma_jox_sinhDec(oxy),...
    'o', 'Color', cs(oxy, :),...
    'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')
end

% do approximate linear fit of lambda_jox(lambda_coxy)

%set param increment
dp = [0.001 0.001];
%set param init guess
p0 = [1 1];
%set param boundaries
pmin = [0.01 0.001];
pmax = [200 2000];
%empty optional params
c = [];

fit_start = 3;
fit_end = 10;
max_iter = 1000;

y_weights = 1./sigma_jox_sinhDec;
[p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                 lm(@linear_model, p0,...
                 coxy_sinhDec', jox_sinhDec', y_weights(fit_start:fit_end), dp,...
                 pmin, pmax, c, fit_start, fit_end, max_iter);
% set fitrange according to values of coxy_sinhDec at fit_start/fit_end

fitrange = linspace(min(coxy_sinhDec), max(coxy_sinhDec), 100);

plot(fitrange, linear_model(fitrange, p), 'LineWidth', 1.5, 'Color','red', ...
    'DisplayName','lin. fit slope m='+string(p(1)))

ylim([5, 30])

xlim([18, 110])
set(gca,'FontSize',15)
xlabel('\lambda_{c} (\mu m)')
ylabel('\lambda_{J_{ox}} (\mu m)')
title('J_{ox} versus c gradient decay length', 'Interpreter','tex')

cb = colorbar;
clim([1 16])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 16, 16);
cb.TickLabels = round(precise_oxygen_levels, 2);
legend('Location','southeast')

savefig(string(plot_path)+'LambdaJox_vs_lambdaCox'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'LambdaJox_vs_lambdaCox'+string(temp_string)+'.png')

%% double-log plot of dependency of lambda_coxy on external oxygen
fig = figure(2);
hold on
errorbar(precise_oxygen_levels, coxy_sinhDec, sigma_coxy_sinhDec, sigma_coxy_sinhDec,...
         'o', 'Color','red', ...
         'MarkerSize',15, 'LineWidth',1.5, 'DisplayName','data')
plot(precise_oxygen_levels, 6*precise_oxygen_levels.^0.6, ...
     'LineWidth',1.5, 'LineStyle','--', 'Color','black', 'DisplayName','~c_{oxy, ext}^{0.6}')

xscale('log')
yscale('log')

legend('Location','southeast')
set(gca,'FontSize',15)
xlabel('c_{oxy, ext.} (\mu M)')
ylabel('\lambda_{c_{oxy}}(c_{oxy, ext})')
title('double log plot of \lambda_{c_{oxy}} versus outside oxygen', 'Interpreter','tex')
savefig(string(plot_path)+'logLambdaCox_vs_logExtOxy'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'logLambdaCox_vs_logExtOxy'+string(temp_string)+'.png')