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



%% plot of ring-averaged (whole-cell) Jox as a function of temperature

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600],...
    'Color','white');
fs = 28;

set(gca,'FontSize',fs)

% get whole-cell Jox by weighted average over Jox(r) for all temperatures
temps = [22 28 31 36];
ring_averaged_jox_list_low = [];
ring_averaged_sigma_jox_list_low = [];

ring_averaged_jox_list_mid = [];
ring_averaged_sigma_jox_list_mid = [];

ring_averaged_jox_list_high = [];
ring_averaged_sigma_jox_list_high = [];

cs = plasma(numel(temps));
colormap(cs)
cs2 = flip(viridis(16));

low_oxy = 5;
mid_oxy = 20;
high_oxy = 50;
offset = 0.8;

for temp=1:numel(temps)
    if temp==1
        vis = 'on';
    else
        vis = 'off';
    end
    [lowv, lowi] = min(abs(precise_oxy_level_list{temp}-low_oxy));
    [midv, midi] = min(abs(precise_oxy_level_list{temp}-mid_oxy));
    [maxv, maxi] = min(abs(precise_oxy_level_list{temp}-high_oxy));

    temp_string = '_T'+string(temps(temp))+'C';
load(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat');
% low oxygen
ring_averaged_jox_list = [];
ring_averaged_sigma_jox_list = [];
for ind=1:(midi-1)
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
ring_averaged_sigma_jox_list_low = [ring_averaged_sigma_jox_list_low mean(ring_averaged_sigma_jox_list)];
ring_averaged_jox_list_low = [ring_averaged_jox_list_low mean(ring_averaged_jox_list)];

% mid oxygen
ring_averaged_jox_list = [];
ring_averaged_sigma_jox_list = [];
for ind=midi:(maxi-1)
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
ring_averaged_sigma_jox_list_mid = [ring_averaged_sigma_jox_list_mid mean(ring_averaged_sigma_jox_list)];
ring_averaged_jox_list_mid = [ring_averaged_jox_list_mid mean(ring_averaged_jox_list)];

% high oxygen
ring_averaged_jox_list = [];
ring_averaged_sigma_jox_list = [];
for ind=maxi:numel(precise_oxygen_levels)
    
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
ring_averaged_sigma_jox_list_high = [ring_averaged_sigma_jox_list_high mean(ring_averaged_sigma_jox_list)];
ring_averaged_jox_list_high = [ring_averaged_jox_list_high mean(ring_averaged_jox_list)];

end

% plotting


hold on

for t=1:4
    if t==1
        vis='on';
    else
        vis='off';
    end
    
% plot(temps(t), ring_averaged_jox_list_low(t), 'o',...
%      'Color', cs2(lowi, :), 'MarkerSize',15, 'LineWidth',1.5,...
%      'DisplayName',"$c_\mathrm{out}\leq$"+string(mid_oxy)+"$\mu$M", 'HandleVisibility',vis)
bar(temps(t), ring_averaged_jox_list_low(t), ...
         'FaceColor', cs2(lowi, :), 'LineWidth',1.5, ...
          'DisplayName',"$c_\mathrm{out}\leq$"+string(mid_oxy)+"$\mu$ M", HandleVisibility=vis)
errorbar(temps(t), ring_averaged_jox_list_low(t), ring_averaged_sigma_jox_list_low(t), ...
         ring_averaged_sigma_jox_list_low(t),...
         'Color', 'black', 'MarkerSize',15, 'LineWidth',1.5, ...
         HandleVisibility='off')


% plot(temps(t)+offset, ring_averaged_jox_list_mid(t), 'o',...
%      'Color', cs2(midi, :), 'MarkerSize',15, 'LineWidth',1.5, ...
%      'DisplayName',string(mid_oxy)+"$\mu$M "+"$\leq c_\mathrm{out}\leq$"+string(high_oxy)+"$\mu$M", 'HandleVisibility',vis)
bar(temps(t)+offset, ring_averaged_jox_list_mid(t),...
         'FaceColor', cs2(midi, :), 'LineWidth',1.5, ...
         'DisplayName',string(mid_oxy)+"$\mu$ M "+"$\leq c_\mathrm{out}\leq$"+string(high_oxy)+"$\mu$ M", 'HandleVisibility',vis)
errorbar(temps(t)+offset, ring_averaged_jox_list_mid(t), ring_averaged_sigma_jox_list_mid(t), ...
         ring_averaged_sigma_jox_list_mid(t),...
         'Color', 'black', 'MarkerSize',15, 'LineWidth',1.5, ...
         'HandleVisibility','off')

% plot(temps(t)+offset*2, ring_averaged_jox_list_high(t), 'o',...
%      'Color', cs2(maxi, :), 'MarkerSize',15, 'LineWidth',1.5, ...
%      'DisplayName',string(high_oxy)+"$\mu$M "+"$\leq c_\mathrm{out}$", 'HandleVisibility',vis)
bar(temps(t)+offset*2, ring_averaged_jox_list_high(t),...
         'FaceColor', cs2(maxi, :), 'LineWidth',1.5, ...
         'HandleVisibility',vis, 'DisplayName',string(high_oxy)+"$\mu$ M "+"$\leq c_\mathrm{out}$")
errorbar(temps(t)+offset*2, ring_averaged_jox_list_high(t), ring_averaged_sigma_jox_list_high(t), ...
         ring_averaged_sigma_jox_list_high(t),...
         'Color', 'black', 'MarkerSize',15, 'LineWidth',1.5, ...
         HandleVisibility='off')
end
xlabel('temperature [$^\circ$C]', 'Interpreter','latex')
%ylabel(' $\overline{\textrm{J}}_{\textrm{ox}}(\textrm{r})$ ($\muM/s$)', 'Interpreter', 'latex')
ylabel('$J_{\textrm{ox}}$ [$\mu$ M/s]', 'Interpreter','latex')
%title('whole-cell average $\mathrm{J}_{\mathrm{ox}}$ versus temperature', 'Interpreter','latex')
legend('Interpreter','latex', 'Location','northwest')
xlim([20 40])
xticks(temperature+offset)
xticklabels(temperature)

savefig(string(plot_path)+'ringAverageJox_vs_temp.fig')
saveas(fig, string(plot_path)+'ringAverageJox_vs_temp.png')
