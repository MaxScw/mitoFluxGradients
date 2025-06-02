%% compare all decay lengths obtained from linear RD, J_ox log-lin plot, c_oxy log-lin plot
clear all
clc

pp_path = '../data/EXP_temperatureResolved/postprocessing_results/';
plot_path = '../data/EXP_temperatureResolved/plots/';

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

    % appen loaded data to list
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
end

%% plot decay length scale of Jox versus that of c_ox for different temps

fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
set(gca,'FontSize',15)
% yscale('log')
ylim([5, 25])
hold on

%cs = viridis(16);
cs = viridis(numel(temps));
colormap(cs)

for t=1:numel(temps)
    for oxy=1:16
        
        if sigma_jox_dec_list{t}(oxy)>3*jox_dec_list{t}(oxy)
            % exception if uncertainty of decay length is very large
            sigma_jox_dec_list{t}(oxy) = 0;
            % errorbar(coxy_dec_list{t}(oxy), jox_dec_list{t}(oxy),...
            % sigma_jox_dec_list{t}(oxy), sigma_jox_dec_list{t}(oxy),...
            % sigma_coxy_dec_list{t}(oxy), sigma_coxy_dec_list{t}(oxy),...
            % 'o', 'Color', 'red',...
            % 'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')
            % errorbar(oxy, jox_dec_list{t}(oxy),...
            % sigma_jox_dec_list{t}(oxy), sigma_jox_dec_list{t}(oxy),...
            % sigma_coxy_dec_list{t}(oxy), sigma_coxy_dec_list{t}(oxy),...
            % 'o', 'Color', 'red',...
            % 'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')

        else
        % plot lambda_jox(lambda_coxy) for different T
        errorbar(coxy_dec_list{t}(oxy), jox_dec_list{t}(oxy),...
            sigma_jox_dec_list{t}(oxy), sigma_jox_dec_list{t}(oxy),...
            sigma_coxy_dec_list{t}(oxy), sigma_coxy_dec_list{t}(oxy),...
            'o', 'Color', cs(t, :),...
            'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')

        % plot lambda_jox(coxy) for different T
        % errorbar(oxy, jox_dec_list{t}(oxy),...
        %     sigma_jox_dec_list{t}(oxy), sigma_jox_dec_list{t}(oxy),...
        %     sigma_coxy_dec_list{t}(oxy), sigma_coxy_dec_list{t}(oxy),...
        %     'o', 'Color', cs(t, :),...
        %     'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')
        end
    end
    
    %(optional) also do approximate linear fit

    % dp = [0.001 0.001];
    % p0 = [1 1];
    % pmin = [0.01 0.001];
    % pmax = [200 2000];
    % c = [];
    
    % % fit
    % 
    % fit_start = 2;
    % fit_end = 9;
    % max_iter = 1000;
    % 
    % y_weights = 1./sigma_linRD_dec;
    % [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
    %                  lm(@linear_model, p0,...
    %                  coxy_sinhDec', linRD_dec', y_weights(fit_start:fit_end), dp,...
    %                  pmin, pmax, c, fit_start, fit_end, max_iter);
    % fitrange = linspace(5, 40, 100);
    % 
    % plot(fitrange, linear_model(fitrange, p), 'LineWidth', 1.5, 'Color','red', ...
    %     'DisplayName','lin. fit slope m='+string(p(1)))
end

set(gca,'FontSize',15)
xlabel('\lambda_{c} (\mu m)')
ylabel('\lambda_{J_{ox}} (\mu m)')
title('J_{ox} versus c gradient decay length', 'Interpreter','tex')

cb = colorbar;
clim([0.5 4.5])
title(cb, 'T (K)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 4, 4);
cb.TickLabels = kelvin_temps;

savefig(string(plot_path)+'LambdaJox_vs_lambdaCox_tempComp.fig')
saveas(fig, string(plot_path)+'LambdaJox_vs_lambdaCox_tempComp.png')

%% plot Jox decay length at maximum oxygen as a function of temperature

fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
set(gca,'FontSize',15)
% yscale('log')
xlim([20+273.5, 40+273.5])
hold on

for temp=1:numel(temps)

    max_jox_dec = jox_dec_list{temp}(end);
    sigma_max_jox_dec = sigma_jox_dec_list{temp}(end);
    plot(kelvin_temps(temp), max_jox_dec, 'o', ...
         'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)
    errorbar(kelvin_temps(temp), max_jox_dec, sigma_max_jox_dec, ...
         'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)

end

set(gca,'FontSize',15)
xlabel('temperature (K)')
ylabel('\lambda_{J_{ox}} (\mu m)')
title('\lambda( J_{ox}) at maximum oxygen versus temperature', 'Interpreter','tex')

savefig(string(plot_path)+'JoxMaxCoxy_vs_temp.fig')
saveas(fig, string(plot_path)+'JoxMaxCoxy_vs_temp.fig')


%%
fig2 = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
hold on
for temp=1:numel(temps)
    % load corrected oxygen levels
    temp_string = '_T'+string(temps(temp))+'C';
    load(string(pp_path)+'cOxy_corr'+string(temp_string)+'.mat')
    
    errorbar(cOxy_all(2, 1:end, end), jox_dec_list{temp}, sigma_jox_dec_list{temp}, ...
        sigma_jox_dec_list{temp}, 'o', 'Color',cs(temp, :), 'MarkerSize',15, 'LineWidth',1.5)
    


end

set(gca,'FontSize',15)
xlabel('external oxygen')
ylabel('\lambda_{J_{ox}} (\mu m)')
%% plot of ring-averaged (whole-cell) Jox as a function of temperature

% get whole-cell Jox by weighted average over Jox(r) for all temperatures
temps = [22 28 31 36];
ring_averaged_jox_list = [];
ring_averaged_sigma_jox_list = [];

for temp_ind=1:numel(temps)
    temp_string = '_T'+string(temps(temp_ind))+'C';
load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
% pick data set for highest external oxygen
ind = 16;
data = [mean(oxygen_ranges_data{ind}.dist_all, 'omitnan');...
        mean(oxygen_ranges_data{ind}.jox_cell_kn_all, 'omitnan')];
data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), 'omitnan');...
               std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), 'omitnan')];

inner_radius = data(1, :);
inner_radius(2:end) = data(1, 1:end-1);
inner_radius(1) = 0;
ring_weights = ((data(1, :).^2 - inner_radius.^2))./(data(1, end).^2);
ring_averaged_jox = sum(ring_weights.*data(2, :));
ring_averaged_sigma_jox = sum(ring_weights.*data_stderr(2, :));

ring_averaged_jox_list = [ring_averaged_jox_list ring_averaged_jox];
ring_averaged_sigma_jox_list = [ring_averaged_sigma_jox_list ring_averaged_sigma_jox];
end

% plotting
fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
set(gca,'FontSize',15)
xlim([273.5+21, 273.5+37])
hold on

plot(kelvin_temps, ring_averaged_jox_list, 'o',...
     'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)
errorbar(kelvin_temps, ring_averaged_jox_list, ring_averaged_sigma_jox_list, ...
         ring_averaged_sigma_jox_list, 'o',...
         'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)

xlabel('temperature (K)')
ylabel('\langle J_{ox}(r) \rangle_r (\muM/s)')
title('whole-cell average J_{ox} versus temperature', 'Interpreter','tex')

savefig(string(plot_path)+'ringAverageJox_vs_temp.fig')
saveas(fig, string(plot_path)+'ringAverageJox_vs_temp.png')

%% plot eff. reaction rate as function of inverse temperature
fig2 = figure('Renderer', 'painters', 'Position', [10 10 1200 600]);
start_ring = 2;
cs = parula(10);
colormap(cs)

load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
dist = mean(oxygen_ranges_data{1}.dist_all, 'omitnan');

km = true;

deltaG_bykb = [];

for ring=start_ring:10
    vMax_temps = [];
    sigma_vMax_temps = [];
    kM_temps = [];
    sigma_kM_temps = [];

    for temp=1:numel(temps)
   
        kM = kM_list{temp}.kM_all(2, :);
        sigma_kM = sigma_kM_list{temp}.sigma_kM_all(2, :);
    
        vMax = vMax_list{temp}.vMax_all(2, :);
        sigma_vMax = sigma_vMax_list{temp}.sigma_vMax_all(2, :);
     
        vMax_temps = [vMax_temps, vMax(ring)];
        sigma_vMax_temps = [sigma_vMax_temps, sigma_vMax(ring)];
        kM_temps = [kM_temps, kM(ring)];
        sigma_kM_temps = [sigma_kM_temps, sigma_kM(ring)];
        
        

    end
    
    subplot(1, 2, 1)
    set(gca,'FontSize',15)
    xlabel('1/T (1/K)')
    
    title('eff. reaction rate per ring')
    %xlim([294, 312])
    hold on
    
    if km == true
    ylabel('log(v_{max}/T)')
    log_rate = log(vMax_temps./(kelvin_temps));
    sigma_log_rate = sqrt((sigma_vMax_temps./vMax_temps).^2 );
    else
    ylabel('log(k_M/T)')
    log_rate = log(kM_temps./(kelvin_temps));
    sigma_log_rate = sqrt((sigma_kM_temps./kM_temps).^2);
    end
    errorbar(1./kelvin_temps, log_rate, sigma_log_rate,...
    'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(ring, :))

    y_weights = 1./sigma_log_rate;

    % fit exp model to v_max(r)
    % param = [slope, intercept]
    
    dp = [0.01 0.001];
    
    p0 = [1e6 2];
    
    pmin = [-1e5 -1e5];
    pmax = [1e5 1e5];
    
    % fit
    
    fit_start = 1;
    fit_end = numel(sigma_log_rate);
    max_iter = 1000;
   
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@linear_model, p0,...
                     100.*1./kelvin_temps', log_rate', y_weights(fit_start:fit_end)', dp,...
                     pmin, pmax, [], fit_start, fit_end, max_iter);
    
    plot(1./kelvin_temps(fit_start:fit_end), linear_model(100.*1./kelvin_temps(fit_start:fit_end), p),...
    'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(ring, :))
    p
    deltaG_bykb = [deltaG_bykb, -p(1)];

    subplot(1, 2, 2)
    xlabel('r (\mu m)')
    ylabel('\Delta G (J/mol)')
    title('eff. activation energy')
    set(gca,'FontSize',15)
    hold on
    errorbar(dist(ring), -p(1)*100*8.3, -sigma_p(1)*100*8.3,...
        'o','MarkerSize',15, 'LineWidth',1.5, 'Color',cs(ring, :))
    %errorbar(1./kelvin_temps, kM_temps, sigma_vMax_temps)
        
end

cb = colorbar;
clim([0.5 10.5-start_ring])
title(cb, 'r (\mu m)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 10-start_ring, 10-start_ring);
cb.TickLabels = dist(start_ring:end);

if km == true
savefig(string(plot_path)+'ringWiseKmReactionRate_vs_temp.fig')
saveas(fig2, string(plot_path)+'ringWiseKmReactionRate_vs_temp.png')
else
savefig(string(plot_path)+'ringWiseVmaxReactionRate_vs_temp.fig')
saveas(fig2, string(plot_path)+'ringWiseVmaxReactionRate_vs_temp.png')
end
%%
p
%% plot vMax and kM profiles for different temperatures

% optional: rescale v_max by 1/v_max(end) and 1/T, kM by 1/T
rescaled_params = false;

start_ring = 2;
fig = figure('Renderer', 'painters', 'Position', [10 10 1000 600]);

set(gca,'FontSize',15)
cs = viridis(numel(temps));
colormap(cs)

hold on

load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
dist = mean(oxygen_ranges_data{1}.dist_all, 'omitnan');

vMaxFit_temp_list = [];
sigma_vMaxFit_temp_list = [];
kMFit_temp_list = [];
sigma_kMFit_temp_list = [];

for temp=1:numel(temps)
    subplot(1, 2, 1)
    set(gca,'FontSize',15)
    hold on
    if rescaled_params==true
        vMax = vMax_list{temp}.vMax_all(2, start_ring:end)./(vMax_list{temp}.vMax_all(2, end)*temps(temp));
        ylabel('v_{max}/(v_{max}(R) T) (\mu M/s)')
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

        ylabel('v_{max} (\mu M/s)')
    end
    errorbar(dist(start_ring:end), vMax, sigma_vMax,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    dist_range = linspace(dist(start_ring), dist(end), 100);
    plot(dist_range, exp_model(dist_range, p),...
        'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    
    xlim([0, 36])
    xlabel('distance (\mu m)')

    subplot(1, 2, 2)
    set(gca,'FontSize',15)
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

        ylabel('k_M (\mu M)')
    end
    errorbar(dist(start_ring:end), kM, sigma_kM,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    dist_range = linspace(dist(start_ring), dist(end), 100);
    plot(dist_range, linear_model(dist_range, p),...
        'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    xlim([0, 36])
    xlabel('distance (\mu m)')
    
end

cb = colorbar;
clim([0.5 4.5])
title(cb, 'T (K)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 4, 4);
cb.TickLabels = kelvin_temps;

if rescaled_params==true
    savefig(string(plot_path)+'vMaxKm_vs_dist_temp_rescaled.fig')
    saveas(fig, string(plot_path)+'vMaxKm_vs_dist_temp_rescaled.png')
else
    savefig(string(plot_path)+'vMaxKm_vs_dist_temp.fig')
    saveas(fig, string(plot_path)+'vMaxKm_vs_dist_temp.png')
end

savefig(string(plot_path)+'JoxmmParamProfile_vs_temp.fig')
saveas(fig, string(plot_path)+'JoxmmParamProfile_vs_temp.png')

%% visualise dependency of k_m(r) fit parameters on temperature

fig = figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
set(gca,'FontSize',15)

for temp_ind=1:numel(temps)
    subplot(1, 2, 1)
    set(gca,'FontSize',15)
    xlim([273.5+20 270+40])
    xlabel("T (K)")
    ylabel("k_M(r) linear fit offset (\mu M)")
    hold on
    temp = kelvin_temps(temp_ind);
    % slope
    errorbar(temp, kMFit_temp_list(1, temp_ind),sigma_kMFit_temp_list(1, temp_ind),...
          'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp_ind, :))
    
    subplot(1, 2, 2)
    set(gca,'FontSize',15)
    xlim([273.5+20 270+40])
    xlabel("T (K)")
    ylabel("k_M(r) linear fit slope (\mu M / \mu m)")
    hold on
    % offset
    errorbar(temp, kMFit_temp_list(2, temp_ind),sigma_kMFit_temp_list(2, temp_ind),...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp_ind, :))
end

%% visualise dependency of v_max fit parameters on temperature

fig = figure('Renderer', 'painters', 'Position', [10 10 1800 600]);
set(gca,'FontSize',15)

for temp_ind=1:numel(temps)
    subplot(1, 3, 1)
    set(gca,'FontSize',15)
    xlim([273.5+20 270+40])
    xlabel("T (K)")
    ylabel("v_{max}(r) exp fit amplitude (\mu M/s)")
    hold on
    temp = kelvin_temps(temp_ind);
    % amplitude
    errorbar(temp, vMaxFit_temp_list(1, temp_ind),sigma_vMaxFit_temp_list(1, temp_ind),...
          'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp_ind, :))
    
    subplot(1, 3, 2)
    set(gca,'FontSize',15)
    xlim([273.5+20 270+40])
    xlabel("T (K)")
    ylabel("v_{max}(r) exp fit decay length (\mu m)")
    hold on
    % decay length
    errorbar(temp, 1./vMaxFit_temp_list(2, temp_ind),sigma_vMaxFit_temp_list(2, temp_ind)./vMaxFit_temp_list(2, temp_ind).^2,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp_ind, :))

    subplot(1, 3, 3)
    set(gca,'FontSize',15)
    xlim([273.5+20 270+40])
    xlabel("T (K)")
    ylabel("v_{max}(r) exp fit offset (\mu M/s)")
    hold on
    % offset
    errorbar(temp, vMaxFit_temp_list(3, temp_ind),sigma_vMaxFit_temp_list(3, temp_ind),...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp_ind, :))
end

corr_vMax_amplitude = vMaxFit_temp_list(3, :);