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





%% plot lambda(J_ox) as a function of lambda(c_oxy) for sinh fit method
fig = figure('Renderer', 'painters', 'Position', [10 10 900 600], 'Color','white');
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
dp = [0.001 0.00];
%set param init guess
p0 = [1 0];
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

fitrange = linspace(0, max(coxy_sinhDec), 100);

plot(fitrange, linear_model(fitrange, [1; 0]), 'LineWidth', 1.5, 'Color','black', ...
    'LineStyle', '--','DisplayName','slope m=1')%+string(round(p(1), 3)))

ylim([0, 110])

xlim([0, 110])
set(gca,'FontSize',19)
xlabel('$\lambda_{\hat{c}}$ [\,\,\,\,m]', 'Interpreter','latex')
ylabel('$\lambda_{J_{\mathrm{ox}}}$ [\,\,\,\,m]', 'Interpreter','latex')
%title('grad. dec. lengths of $J_{\mathrm{ox}}$ and $\hat{c}$', 'Interpreter','latex')

% cb = colorbar;
% clim([1 16])
% title(cb, 'c_\mathrm{out} (\mu M)', 'Interpreter', 'tex')
% 
% clim([1.5 numel(precise_oxygen_levels)+0.5])
% title(cb, '$c_\mathrm{out} [\mu M]$', 'Interpreter', 'latex')
% cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
% cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);

% cb.Ticks = linspace(1, 16, 16);
% cb.TickLabels = round(precise_oxygen_levels, 2);
legend('Location','southeast', 'Interpreter','latex')

savefig(string(plot_path)+'LambdaJox_vs_lambdaCox'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'LambdaJox_vs_lambdaCox'+string(temp_string)+'.png')
export_fig(fig, string(plot_path)+'LambdaJox_vs_lambdaCox'+string(temp_string)+'.eps')