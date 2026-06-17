%% load setup file for file paths and data
clc
clear all
temp_resolved = false;
run("setup.m")

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
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],...
    'Color','white');
fs = 19;

oxy_range = linspace(0, oxygen_levels(end), 100);
colors = flip(cool(15));
colormap(colors(1:10, 1:end));

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
errorbar(corr_oxy(1:end, 2), mean(flux, 1),...
         mean(sigma_flux, 1), mean(sigma_flux, 1), 'o', 'Color', colors(2, :), ...
         'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off')
ax=gca;
set(gca,'FontSize',fs);

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
errorbar(corr_oxy(1:end, 8), mean(flux, 1),...
         mean(sigma_flux, 1), mean(sigma_flux, 1), 'o', 'Color', colors(8, :), ...
         'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off')
hold on;
ax=gca;
set(gca,'FontSize',fs);

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
    'DisplayName','$J_{\mathrm{ox}}^{\mathrm{fit}}$')
errorbar(corr_oxy(1:end, 10), mean(flux, 1),...
         mean(sigma_flux, 1), mean(sigma_flux, 1), 'o', 'Color', colors(10, :), ...
         'MarkerSize',10, 'LineWidth',1.5, 'DisplayName', '$J_{\mathrm{ox}}$')
hold on;
ax=gca;
set(gca,'FontSize',fs);
legend('Interpreter','latex')
xlabel('$\hat{c}(c_\mathrm{out}, r)$ [\,\,\,\,M]', 'Interpreter','latex')
ylabel('ETC flux [\,\,\,\,M/s]', 'Interpreter','latex')
%title('michaelis-menten rate law fit', 'Interpreter','latex')
xlim([-5 55])
ylim([-4 135])
cb = colorbar;
clim([0.5 10.5])
title(cb, '$r$ [\,\,\,\,m]', 'Interpreter', 'latex')

cb.Ticks = linspace(1, 10, 10);
cb.TickLabels = round(mean(dist_per_ring, 2), 2);

savefig(string(plot_path)+'Jox_vs_coxy_mmFit'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_coxy_mmFit'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'Jox_vs_coxy_mmFit'+string(temp_string)+'.eps')