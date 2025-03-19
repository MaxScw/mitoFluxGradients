%% compare all decay lengths obtained from linear RD, J_ox log-lin plot, c_oxy log-lin plot
clear all
clc

pp_path = '../data/EXP_temperatureResolved/postprocessing_results/';
plot_path = '../data/EXP_temperatureResolved/plots/';

temps = [22 28 31 36];

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
    
    % load J_ox decay lengths from log lin fit
    % joxEmp = load('decayLength_joxEmp_T'+string(temps(t))+'C.mat');
    % sigma_joxEmp = load('sigma_decayLength_joxEmp_T'+string(temps(t))+'C.mat');
    % joxEmp_dec = joxEmp.('decay_lengths');
    % sigma_joxEmp_dec = sigma_joxEmp.('sigma_decay_lengths');
    
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
    jox_dec_list{end+1} = linRD_dec;
    sigma_jox_dec_list{end+1} = sigma_linRD_dec;
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
title(cb, 'T (°C)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 4, 4);
cb.TickLabels = temps;

savefig(string(plot_path)+'LambdaJox_vs_lambdaCox_tempComp.fig')
saveas(fig, string(plot_path)+'LambdaJox_vs_lambdaCox_tempComp.fig')

%% plot Jox decay length at maximum oxygen as a function of temperature

fig = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
set(gca,'FontSize',15)
% yscale('log')
xlim([20, 40])
hold on

for temp=1:numel(temps)

    max_jox_dec = jox_dec_list{temp}(end);
    sigma_max_jox_dec = sigma_jox_dec_list{temp}(end);
    plot(temps(temp), max_jox_dec, 'o', ...
         'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)
    errorbar(temps(temp), max_jox_dec, sigma_max_jox_dec, ...
         'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)

end

set(gca,'FontSize',15)
xlabel('temperature (°C)')
ylabel('\lambda_{J_{ox}} (\mu m)')
title('\lambda( J_{ox}) at maximum oxygen versus temperature', 'Interpreter','tex')

savefig(string(plot_path)+'JoxMaxCoxy_vs_temp.fig')
saveas(fig, string(plot_path)+'JoxMaxCoxy_vs_temp.fig')

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
xlim([21, 37])
hold on

plot(temps, ring_averaged_jox_list, 'o',...
     'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)
errorbar(temps, ring_averaged_jox_list, ring_averaged_sigma_jox_list, ...
         ring_averaged_sigma_jox_list, 'o',...
         'Color', 'red', 'MarkerSize',15, 'LineWidth',1.5)

xlabel('temperature (°C)')
ylabel('\langle J_{ox}(r) \rangle_r (\muM/s)')
title('whole-cell average J_{ox} versus temperature', 'Interpreter','tex')

savefig(string(plot_path)+'ringAverageJox_vs_temp.fig')
saveas(fig, string(plot_path)+'ringAverageJox_vs_temp.png')

%% plot vMax and kM profiles for different temperatures

% optional: rescale v_max by 1/v_max(end) and 1/T, kM by 1/T
rescaled_params = true;

start_ring = 2;
fig = figure('Renderer', 'painters', 'Position', [10 10 1000 600]);

set(gca,'FontSize',15)
cs = viridis(numel(temps));
colormap(cs)

hold on

load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
dist = mean(oxygen_ranges_data{1}.dist_all, 'omitnan');

for temp=1:numel(temps)
    subplot(1, 2, 1)
    set(gca,'FontSize',15)
    hold on
    if rescaled_params==true
        vMax = vMax_list{temp}.vMax_all(2, start_ring:end)./(vMax_list{temp}.vMax_all(2, end)*temps(temp));
        ylabel('v_{max}/(v_{max}(R) T) (\mu M/s)')
    else
        vMax = vMax_list{temp}.vMax_all(2, start_ring:end);
        ylabel('v_{max} (\mu M/s)')
    end
    plot(dist(start_ring:end), vMax,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
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
        ylabel('k_M (\mu M)')
    end
    plot(dist(start_ring:end), kM,...
        'o', 'MarkerSize',15, 'LineWidth',1.5, 'Color',cs(temp, :))
    xlim([0, 36])
    xlabel('distance (\mu m)')
    
end

cb = colorbar;
clim([0.5 4.5])
title(cb, 'T (°C)', 'Interpreter', 'tex')

cb.Ticks = linspace(1, 4, 4);
cb.TickLabels = temps;

if rescaled_params==true
    savefig(string(plot_path)+'vMaxKm_vs_dist_temp_rescaled.fig')
    saveas(fig, string(plot_path)+'vMaxKm_vs_dist_temp_rescaled.png')
else
    savefig(string(plot_path)+'vMaxKm_vs_dist_temp.fig')
    saveas(fig, string(plot_path)+'vMaxKm_vs_dist_temp.png')
end    