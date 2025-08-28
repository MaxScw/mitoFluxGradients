%% designation of oxygen ranges to analyse
clc
clear all

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

%% loading of mito density, temperatures, solubilities, diffusivities

load(data_path + "mito_density_oocytes/mito_density.mat")
mean_mito_density = mean(ratio_num_all_mitotracker, 1);

load(data_path + 'exp_solubilities.mat')
load(data_path + 'exp_temperatures.mat')
load(data_path + "exp_diffusivities.mat")

if temp_resolved==true
    temp_ind = 2;
    temp_string = '_T'+string(temperature(temp_ind))+'C';
else 
    temp_string = '';
end

%% loading of read-in data

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
end
%% load corrected internal oxygen
corr_oxy = squeeze(load(string(pp_path)+'cOxy_corr'+string(temp_string)+'.mat').cOxy_all(2, :, :));

%% fit michealis menten rate law to J_ox(c_oxy) per ring
jox_per_ring = [];
sigma_jox_per_ring = [];
dist_per_ring = [];

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];

for ind=1:numel(oxygen_ranges)

    data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];
    jox_per_ring = [jox_per_ring, data(2, :)'];
    sigma_jox_per_ring = [sigma_jox_per_ring, data_stderr(2, :)'];
    dist_per_ring = [dist_per_ring, data(1, :)'];
end

for ring=1:10
    jox = jox_per_ring(ring, :);
    sigma_jox = sigma_jox_per_ring(ring, :);
    y_weights = 1./sigma_jox;

    oxygen_levels = corr_oxy(:, ring);
    
    % set parameters to fit, increment and initial values (LAMBDA LINEAR MODEL)
    % param = [vmax, km]
    
    dp = [0.001 0.001];

    p0 = [jox(end) precise_oxygen_levels(end)];

    pmin = [0 1e-5];
    pmax = [200 200];
    
    % fit
    
    fit_start = 2;
    fit_end = size(precise_oxygen_levels, 2);
    max_iter = 10000;
    
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@mm_flux, p0,...
                     oxygen_levels, jox', y_weights(fit_start:fit_end), dp,...
                     pmin, pmax, [], fit_start, fit_end, max_iter);
    fit_params = [fit_params, p];
    sigma_fit_params = [sigma_fit_params, sigma_p];
    chi_squareds = [chi_squareds X2];
    R_squareds = [R_squareds R_sq];

end

%% plot michaelis-menten rate law for all rings
fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);

oxy_range = linspace(0, oxygen_levels(end), 100);
colors = flip(cool(15));
colormap(colors(1:10, 1:end))

% plot average of J_ox(c_oxy(r)) michaelis-menten fits for three different
% distance regimes (low, mid, high)

%low
flux = [];
sigma_flux = [];
fitflux = [];
hold on;
%average
for ring=2
    flux = [flux; jox_per_ring(ring, :)];
    sigma_flux = [sigma_flux; sigma_jox_per_ring(ring, :)];
    fitflux = [fitflux; mm_flux(oxy_range, fit_params(:, ring), [])];
   
end
plot(oxy_range, mean(fitflux, 1), 'Color', colors(2, :), 'LineWidth',1.5, ...
    'HandleVisibility','off')
errorbar(mean(corr_oxy, 2), mean(flux, 1),...
         mean(sigma_flux, 1), mean(sigma_flux, 1), 'o', 'Color', colors(2, :), ...
         'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off')
ax=gca;
set(gca,'FontSize',19);

%mid
flux = [];
sigma_flux = [];
fitflux = [];
hold on
%average
for ring=8
    flux = [flux; jox_per_ring(ring, :)];
    sigma_flux = [sigma_flux; sigma_jox_per_ring(ring, :)];
    fitflux = [fitflux; mm_flux(oxy_range, fit_params(:, ring), [])];

end
plot(oxy_range, mean(fitflux, 1), 'Color', colors(8, :), 'LineWidth',1.5, ...
    'HandleVisibility','off')
errorbar(mean(corr_oxy, 2), mean(flux, 1),...
         mean(sigma_flux, 1), mean(sigma_flux, 1), 'o', 'Color', colors(8, :), ...
         'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off')
hold on;
ax=gca;
set(gca,'FontSize',19);

%high
flux = [];
sigma_flux = [];
fitflux = [];
%average
for ring=10
    flux = [flux; jox_per_ring(ring, :)];
    sigma_flux = [sigma_flux; sigma_jox_per_ring(ring, :)];
    fitflux = [fitflux; mm_flux(oxy_range, fit_params(:, ring), [])];
end
plot(oxy_range, mean(fitflux, 1), 'Color', colors(10, :), 'LineWidth',1.5, ...
    'DisplayName','J_{ox}^{theory}')
errorbar(mean(corr_oxy, 2), mean(flux, 1),...
         mean(sigma_flux, 1), mean(sigma_flux, 1), 'o', 'Color', colors(10, :), ...
         'MarkerSize',10, 'LineWidth',1.5, 'DisplayName', 'J_{ox}')
hold on;
ax=gca;
set(gca,'FontSize',19);
legend()
xlabel('$\hat{c}(r) (\mu M)$', 'Interpreter','latex')
ylabel('ETC flux $(\mu M/s)$', 'Interpreter','latex')
title('michaelis-menten rate law fit', 'Interpreter','latex')
xlim([-5 55])
ylim([-4 135])
cb = colorbar;
clim([0.5 10.5])
title(cb, 'r $(\mu m)$', 'Interpreter', 'latex')
dist_per_ring(:, 1)
cb.Ticks = linspace(1, 10, 10);
cb.TickLabels = round(mean(dist_per_ring, 2), 2);

savefig(string(plot_path)+'Jox_vs_coxy_mmFit'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_coxy_mmFit'+string(temp_string)+'.png')

%% plot coeff. of determination from mm fit for different oxygen levels

fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
% quality of fit as function of distance
dist = mean(oxygen_ranges_data{1}.dist_all, 'omitnan');
plot(dist(2:end), chi_squareds(2:end), 'o', 'LineWidth',1.5, 'Color','red')
xlabel('distance from cell center (\mu m)', 'Interpreter','tex')
ylabel('R^2 coeff. of determination', 'Interpreter','tex')
ax=gca;
set(gca,'FontSize',15);
title('quality of fit')

savefig(string(plot_path)+'Jox_vs_coxy_mmFitQual'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_coxy_mmFitQual'+string(temp_string)+'.png')