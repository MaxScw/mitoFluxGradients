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

%% plot Jox decay length at maximum oxygen as a function of temperature

fig = figure('Renderer', 'painters', 'Position', [10 10 600 600],...
    'Color','white');
fs = 28;
set(gca,'FontSize',fs)
% yscale('log')
%xlim([20+273.5, 40+273.5])
ylim([5, 30])
hold on

cs = magma(numel(temps));
colormap(cs)
cs2 = flip(viridis(16));

low_oxy = 5;
mid_oxy = 20;
high_oxy = 50;
offset = 0.4;

for temp=1:numel(temps)
    if temp==1
        vis = 'on';
    else
        vis = 'off';
    end
    [lowv, lowi] = min(abs(precise_oxy_level_list{temp}-low_oxy));
    [midv, midi] = min(abs(precise_oxy_level_list{temp}-mid_oxy));
    [maxv, maxi] = min(abs(precise_oxy_level_list{temp}-high_oxy));

    min_jox_dec = mean(jox_dec_list{temp}(1:midi-1));
    mid_jox_dec = mean(jox_dec_list{temp}(midi:maxi-1));
    max_jox_dec = mean(jox_dec_list{temp}(maxi:end));
    sigma_min_jox_dec = mean(sigma_jox_dec_list{temp}(1:midi-1));
    sigma_mid_jox_dec = mean(sigma_jox_dec_list{temp}(midi:maxi-1));
    sigma_max_jox_dec = mean(sigma_jox_dec_list{temp}(maxi:end));
    
    % plot(temps(temp)+offset*2, min_jox_dec, 'o', ...
    %      'Color', cs2(lowi, :), 'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')
    bar(temps(temp)+offset*(-2), min_jox_dec, ...
         'FaceColor', cs2(lowi, :), 'LineWidth',1.5, ...
         'DisplayName',"$c_\mathrm{out}\leq$"+string(mid_oxy)+"$\mu$ M", 'HandleVisibility',vis)
    errorbar(temps(temp)+offset*(-2), min_jox_dec, sigma_max_jox_dec, ...
         'Color', 'black', 'MarkerSize',15, 'LineWidth',1.5, ...
         'DisplayName',"$c_\mathrm{out}\leq$"+string(mid_oxy)+"$\mu$ M", 'HandleVisibility','off')
    
    % plot(temps(temp), mid_jox_dec, 'o', ...
    %      'Color', cs2(midi, :), 'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')
    bar(temps(temp), mid_jox_dec,...
         'FaceColor', cs2(midi, :), 'LineWidth',1.5, ...
         'DisplayName',string(mid_oxy)+"$\mu$ M "+"$\leq c_\mathrm{out}\leq$"+string(high_oxy)+"$\mu$ M", ...
         'HandleVisibility',vis)
    errorbar(temps(temp), mid_jox_dec, sigma_mid_jox_dec, ...
         'Color', 'black', 'MarkerSize',15, 'LineWidth',1.5, ...
         'DisplayName',string(mid_oxy)+"$\mu$ M "+"$\leq c_\mathrm{out}\leq$"+string(high_oxy)+"$\mu$ M", ...
         'HandleVisibility','off')

    % plot(temps(temp)+offset*(-2), max_jox_dec, 'o', ...
    %      'Color', cs2(maxi, :), 'MarkerSize',15, 'LineWidth',1.5, 'HandleVisibility','off')
    bar(temps(temp)+offset*(2), max_jox_dec, ...
     'FaceColor', cs2(maxi, :), 'LineWidth',1.5, ...
     'DisplayName',string(high_oxy)+"$\mu$ M "+"$\leq c_\mathrm{out}$", 'HandleVisibility',vis)
    errorbar(temps(temp)+offset*(2), max_jox_dec, sigma_max_jox_dec, ...
         'Color', 'black', 'MarkerSize',15, 'LineWidth',1.5, ...
         'DisplayName',string(high_oxy)+"$\mu$ M "+"$\leq c_\mathrm{out}$", 'HandleVisibility','off')

end
legend('Location','southeast', 'Interpreter','latex')
set(gca,'FontSize',fs)
xlabel('temperature [$^\circ$C]', 'Interpreter','latex')
ylabel('$\lambda_{J_{\mathrm{ox}}}$ [$\mu$ m]', 'Interpreter','latex')
%title('$\lambda( J_{\mathrm{ox}})$ versus temperature', 'Interpreter','latex')
ylim([0 30])
xticks(temps)
xticklabels(temps)

savefig(string(plot_path)+'Jox_vs_temp.fig')
saveas(fig, string(plot_path)+'Jox_vs_temp.png')



