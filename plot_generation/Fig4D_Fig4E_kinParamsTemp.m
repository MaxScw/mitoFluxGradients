%% load setup file for file paths and data
clc
clear all
temp_resolved = true;
run("setup.m")

%% compare all decay lengths obtained from linear RD, J_ox log-lin plot, c_oxy log-lin plot

temps = [22 28 31 36];
kelvin_temps = temps+273.15;

jox_dec_list = {};
sigma_jox_dec_list = {};
coxy_dec_list = {};
sigma_coxy_dec_list = {};
vMax_list = {};
sigma_vMax_list = {};
kM_list = {};
sigma_kM_list = {};
precise_oxy_level_list = {};


for t=1:numel(temps)
    temp_string = '_T'+string(temps(t))+'C';

    % load J_ox decay lengths from sinh fit
    linRD_fit = load(string(pp_path)+'gradient_fit_linear_model'+string(temp_string)+'.mat');
    linRD_fit = linRD_fit.('fit_data_gradient');
    linRD_fit_params = linRD_fit{1};
    linRD_fit_sigma_params = linRD_fit{2};

    linRD_dec = linRD_fit_params(3, :);
    sigma_linRD_dec = linRD_fit_sigma_params(3, :);
    
    % load J_ox decay length obtained by sinh fit
    jox_fit = load(string(pp_path)+'sinhDecayLength_jox'+string(temp_string)+'.mat');
    jox_sigmafit = load(string(pp_path)+'sigma_sinhDecayLength_jox'+string(temp_string)+'.mat');
    jox_sinhDec = jox_fit.("decay_lengths");
    sigma_jox_sinhDec = jox_sigmafit.("sigma_decay_lengths");
    
    % load c_oxy decay lengths from sinh fit
    coxy_sinhFit = load(string(pp_path)+'sinhDecayLength_coxy'+string(temp_string)+'.mat');
    coxy_sigmaSinhFit = load(string(pp_path)+'sigma_sinhDecayLength_coxy'+string(temp_string)+'.mat');
    coxy_sinhDec = coxy_sinhFit.('decay_lengths');
    sigma_coxy_sinhDec = coxy_sigmaSinhFit.('sigma_decay_lengths');

    % load vMax and kM parameters from mm fit
    kM = load(string(pp_path)+'k_m_profiles_corrected'+string(temp_string)+'.mat');
    sigma_kM = load(string(pp_path)+'k_m_profiles_corrected_sigma'+string(temp_string)+'.mat');
    vMax = load(string(pp_path)+'v_max_profiles_corrected'+string(temp_string)+'.mat');
    sigma_vMax = load(string(pp_path)+'v_max_profiles_corrected_sigma'+string(temp_string)+'.mat');
    


    load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
    
    precise_oxygen_levels = [];
    sigma_precise_oxygen_levels = [];
    for oxy_ind=1:16
        precise_oxygen_levels = [precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.o2_levels];
        sigma_precise_oxygen_levels = [sigma_precise_oxygen_levels, oxygen_ranges_data{oxy_ind}.sigma_o2_levels];
    end
    
 
    precise_oxygen_levels = solubility(t).*precise_oxygen_levels./20.946;
    sigma_precise_oxygen_levels = solubility(t).*sigma_precise_oxygen_levels./20.946;
    sigma_precise_oxygen_levels(sigma_precise_oxygen_levels==0) = mean(sigma_precise_oxygen_levels);



    % append loaded data to list
    jox_dec_list{end+1} = jox_sinhDec;
    sigma_jox_dec_list{end+1} = sigma_jox_sinhDec;
    % jox_dec_list{end+1} = joxEmp_dec;
    % sigma_jox_dec_list{end+1} = sigma_joxEmp_dec;

    coxy_dec_list{end+1} = coxy_sinhDec;
    sigma_coxy_dec_list{end+1} = sigma_coxy_sinhDec;

    vMax_list{end+1} = vMax;
    sigma_vMax_list{end+1} = sigma_vMax;
    kM_list{end+1} = kM;
    sigma_kM_list{end+1} = sigma_kM;

    precise_oxy_level_list{end+1} = precise_oxygen_levels;
end


%% plot vMax profiles for different temperatures

% optional: rescale v_max by 1/v_max(end) and 1/T, kM by 1/T
rescaled_params = false;

load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
dist = mean(oxygen_ranges_data{1}.dist_all, 'omitnan');

start_ring = 2;
fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
    'Color','white');


set(gca,'FontSize',19)
cs = plasma(numel(temps)+1);
cs = cs(1:end-1, 1:end);
colormap(cs)

hold on

vMaxFit_temp_list = [];
sigma_vMaxFit_temp_list = [];
kMFit_temp_list = [];
sigma_kMFit_temp_list = [];

for temp=1:numel(temps)
    
    set(gca,'FontSize',19)
    hold on
    if rescaled_params==true
        vMax = vMax_list{temp}.vMax_all(2, start_ring:end)./(vMax_list{temp}.vMax_all(2, end)*temps(temp));
        ylabel('$v_{\mathrm{max}}/(v_{\mathrm{max}}(\mathrm{R}) \mathrm{T})$ [$\mu$ M/s]', 'Interpreter','latex')
    else
        vMax = vMax_list{temp}.vMax_all(2, start_ring:end);
        sigma_vMax = sigma_vMax_list{temp}.sigma_vMax_all(2, start_ring:end);

        y_weights = 1./sigma_vMax;

        % fit exp model to v_max(r)
        % param = [amp, dec, offset]
        
        dp = [0 0.001 0.001];
        
        p0 = [0.1 1 40];
        
        pmin = [-200 0.001 -200];
        pmax = [200 200 200];
        
        % fit
        
        fit_start = 1;
        fit_end = numel(sigma_vMax);
        max_iter = 1000;
        size(dist(start_ring:end))
        size(vMax)
        size(y_weights)
        
        [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                         lm(@exp_model, p0,...
                         dist(start_ring:end)', vMax', y_weights', dp,...
                         pmin, pmax, [], fit_start, fit_end, max_iter);
        vMaxFit_temp_list = [vMaxFit_temp_list, p];
        sigma_vMaxFit_temp_list = [sigma_vMaxFit_temp_list, sigma_p];
        
        ylabel('$v_{\mathrm{max}}$ [$\mu$ M/s]', 'Interpreter','latex')
    end
    errorbar(dist(start_ring:end), vMax, sigma_vMax,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    dist_range = linspace(dist(start_ring), dist(end), 100);
    % plot(dist_range, exp_model(dist_range, p),...
    %     'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    
    xlim([0, 36])
    xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex')
end

cb = colorbar;
clim([0.5 4.5])
title(cb, 'T [$^\circ$C]', 'Interpreter', 'latex')

cb.Ticks = linspace(1, 4, 4);
cb.TickLabels = temps;

if rescaled_params==true
    savefig(string(plot_path)+'vMax_vs_dist_temp_rescaled.fig')
    saveas(fig, string(plot_path)+'vMax_vs_dist_temp_rescaled.png')
else
    savefig(string(plot_path)+'vMax_vs_dist_temp.fig')
    saveas(fig, string(plot_path)+'vMax_vs_dist_temp.png')
end


fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
    'Color','white');
cs = plasma(numel(temps)+1);
cs = cs(1:end-1, 1:end);
colormap(cs)
set(gca,'FontSize',19)

for temp=1:numel(temps)
    
    
    hold on
    if rescaled_params==true
        kM = kM_list{temp}.kM_all(2, start_ring:end)./temps(temp);
        ylabel('k_M/T (\mu M)')
    else
        kM = kM_list{temp}.kM_all(2, start_ring:end);
        sigma_kM = sigma_kM_list{temp}.sigma_kM_all(2, start_ring:end);

        y_weights = 1./sigma_kM;

        % fit exp model to v_max(r)
        % param = [amp, dec, offset]
        
        dp = [0.001 0.001];
        
        p0 = [10 1/30];
        
        pmin = [-200 -200];
        pmax = [200 200];
        
        % fit
        
        fit_start = 1;
        fit_end = numel(sigma_kM);
        max_iter = 1000;
        [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                         lm(@linear_model, p0,...
                         dist(start_ring:end)', kM', y_weights', dp,...
                         pmin, pmax, [], fit_start, fit_end, max_iter);
        kMFit_temp_list = [kMFit_temp_list, p];
        sigma_kMFit_temp_list = [sigma_kMFit_temp_list, sigma_p];

        ylabel('$K_\mathrm{M}$ [$\mu$ M]', 'Interpreter','latex')
    end
    errorbar(dist(start_ring:end), kM, sigma_kM,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    dist_range = linspace(dist(start_ring), dist(end), 100);
    % plot(dist_range, linear_model(dist_range, p),...
    %     'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    xlim([0, 36])
    xlabel('distance from cell center ($r$) [$\mu$ m]', 'Interpreter','latex')
end



if rescaled_params==true
    savefig(string(plot_path)+'Km_vs_dist_temp_rescaled.fig')
    saveas(fig2, string(plot_path)+'Km_vs_dist_temp_rescaled.png')
else
    savefig(string(plot_path)+'Km_vs_dist_temp.fig')
    saveas(fig2, string(plot_path)+'Km_vs_dist_temp.png')
end
