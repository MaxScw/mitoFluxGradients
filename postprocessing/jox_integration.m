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

% temp_resolved = false;
% data_path = "/home/mx/mitoFluxGradients/data/";
% pp_path = "/home/mx/mitoFluxGradients/data/EXP_published/postprocessing_results/";
% plot_path = "/home/mx/mitoFluxGradients/data/EXP_published/plots/";
temp_resolved = true;
data_path = "/home/mx/mitoFluxGradients/data/";
pp_path = "../data/EXP_temperatureResolved/postprocessing_results/";
plot_path = "../data/EXP_temperatureResolved/plots/";

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

%%
n_rings = 10;

% determine range of diffusion coefficients, for which to calculate

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
    data = [mean(oxygen_ranges_data{ind}.dist_all, "omitnan");...
            mean(oxygen_ranges_data{ind}.jox_cell_kn_all, "omitnan")];
    data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), "omitnan");...
                   std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), "omitnan")];
    jox_data(ind, :) = data(2, :);
    r_data(ind, :) = data(1, :);
    sigma_jox_data(ind, :) = data_stderr(2, :);
    sigma_r_data(ind, :) = data_stderr(1, :);

    % WHETHER MITO DENSITY WEIGHTED JOX DATA IS USED:
    % APPLY mito_weight TO data

    % include inheterogeneity
    mito_weight = mean_mito_density;
    
    % average over inhet.
    %mito_weight = mean(mean_mito_density);

    jox_data_mitoDist(ind, :) = data(2, :).*mito_weight./2;    

end

% output containers for c_oxy, v_max, k_m have dimensions
% 2x(n_oxy_levels)x(n_rings), the first dimension stands for results
% at uncorrected versus corrected intracellular oxygen

cOxy_all = zeros(2, numel(precise_oxygen_levels), n_rings);
cOxy_all(1, :, :) = repmat(precise_oxygen_levels, n_rings, 1)';
kM_all =  zeros(2, n_rings);
sigma_kM_all =  zeros(2, n_rings);
vMax_all =  zeros(2, n_rings);
sigma_vMax_all =  zeros(2, n_rings);
jox_pred = zeros(2, numel(precise_oxygen_levels), n_rings);
%% test: plot a selection of j_ox gradients to check whether mitoch. density weights are applied or not

cs = viridis(16);

for ind=[1 2 3 4]

    % plot(r_data(ind, :), jox_data(ind, :), "o","Color","blue",...
    %      "MarkerSize",10, "LineWidth",1.5)
    hold on
    plot(r_data(ind, :), jox_data_mitoDist(ind, :), "o", "Color",cs(ind, :),...
         "MarkerSize",10, "LineWidth",1.5)
    rrange = linspace(r_data(ind, 1), r_data(ind, end), 100);
    
    plot(rrange, interp1(r_data(ind, :), jox_data_mitoDist(ind, :), rrange), "Color",cs(ind, :),...
         "MarkerSize",10, "LineWidth",1.5)

end

%%
% option to display plot of michaelis-menten fits
print_plot = false;

% option to only record quality of fit for corrected iteration
save_all_its = true;

% create lists for containing quality of fit at certain diff_coeff
x2_hist = [];
r2_hist = [];
r2_jox_hist = [];
x2_jox_hist = [];

% assign first ring index considered for mm param fit & oxy integration
% (for now, v_max, k_m, J_ox and c_oxy values for rings not considered will
% be set to zero)
first_ring = 1;

for diff_ind=1:numel(diff_coeffs)

    for it=1:2
    
        % based on the current c_oxy(r, c*), v_max and k_m are calculated
        % from a mm-fit to J_ox(c_oxy)
        % for iteration 1, c_oxy(r, c*) = c* is assumed
        
        x2_current = [];
        r2_current = [];
        for ring=first_ring:n_rings
            jox = jox_data(:, ring);
            sigma_jox = sigma_jox_data(:, ring);

            % jox = jox_data_mitoDist(:, ring);
            % sigma_jox = jox.*sqrt((sigma_jox_data(:, ring)./jox_data(:, ring)).^2 +...
            %                  (stderr_mito_density./mito_weight).^2);
            c_oxy = cOxy_all(it, :, ring);
    
            %  for each ring, obtain v_max and k_m from mm rate law fit
    
            y_weights = 1./sigma_jox;
            % param = [vmax, km]
            dp = [0.001 0.001];
            p0 = [jox(end) 0.1];
            pmin = [0.001 0.0001];
            pmax = [400 200];
            % fit
            fit_start = 1;
            fit_end = size(precise_oxygen_levels, 2);
            max_iter = 1000;
            squeeze(c_oxy)'
            [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                             lm(@mm_flux, p0,...
                             squeeze(c_oxy)', jox, y_weights(fit_start:fit_end), dp,...
                             pmin, pmax, [], fit_start, fit_end, max_iter);
            v_max = p(1);
            k_m = p(2);
            sigma_v_max = sigma_p(1);
            sigma_k_m = sigma_p(2);
            
            if (it==2)&&(print_plot==true)
            figure(99)
            cs = viridis(n_rings);
            hold on
            % subplot(1, 2, 1)
            % plot((k_m./v_max)*1./c_oxy, 1./jox, "o", "Color",cs(ring, :))
            % xlim([0, 0.1])
            % ylim([0, 0.1])
            % 
            % subplot(1, 2, 2)
            c_oxy_ext = cOxy_all(1, :, ring);
            plot(c_oxy_ext, jox./c_oxy, "o", "Color",cs(ring, :))
            %yscale("log")
            
            %plot(c_oxy, mm_flux(c_oxy, p, []), "Color",cs(ring, :))

            %plot(c_oxy, mm_flux(c_oxy, [14; 0.1], []), "Color",cs(ring, :))
            end

            vMax_all(it, ring) = v_max;
            sigma_vMax_all(it, ring) = sigma_v_max;
            kM_all(it, ring) = k_m;
            sigma_kM_all(it, ring) = sigma_k_m;
    
            x2_current = [x2_current, X2];
            r2_current = [r2_current, R_sq];
        end
        
        if (it==2)|(save_all_its==true)
            x2_hist = [x2_hist, mean(x2_current)];
            r2_hist = [r2_hist, mean(r2_current)];
        end

        % calculate prediction for J_ox(r, c*) based on mm model
        % with varying v_max, k_m (take most recent values)
        r2_jox = [];
        x2_jox = [];
        for ind=1:numel(precise_oxygen_levels)
            
            jox_pred(it, ind, :) = mm_flux(squeeze(cOxy_all(it, ind, :))', [vMax_all(it, :); kM_all(it, :)], []);

            jox_curr = jox_data(ind, :);
            
            R_sq = corrcoef(jox_curr(2:end), squeeze(jox_pred(it, ind, 2:end))');
            R_sq = R_sq(1,2).^2;
            r2_jox = [r2_jox, R_sq];
            
            x2_jox = [x2_jox, sum(((jox_curr(2:end) - squeeze(jox_pred(it, ind, 2:end))').^2)./jox_curr(2:end))];
     
        end
        if (it==2)|(save_all_its==true)
            r2_jox_hist = [r2_jox_hist; r2_jox];
            x2_jox_hist = [x2_jox_hist; x2_jox];
        end
    
        % calculate new oxygen levels by integration of flux equation
        for ind=1:numel(precise_oxygen_levels)
            data = [mean(oxygen_ranges_data{ind}.dist_all, "omitnan");...
                    mean(oxygen_ranges_data{ind}.jox_cell_kn_all, "omitnan")];
            data_stderr = [std(oxygen_ranges_data{ind}.dist_all./sqrt(size(oxygen_ranges_data{ind}.dist_all,1)), "omitnan");...
                           std(oxygen_ranges_data{ind}.jox_cell_kn_all./sqrt(size(oxygen_ranges_data{ind}.jox_cell_kn_all,1)), "omitnan")];
            
            D = diff_coeffs(diff_ind);
            c_star = precise_oxygen_levels(ind);
            params = [c_star, D];
            
            c_oxy_res = flux_integrator(data(1, first_ring:end), params, jox_data_mitoDist(ind, first_ring:end));
            cOxy_all(2, ind, first_ring:end) = c_oxy_res(1, :);
        end
    
    end

end

%%
cs = viridis(16)
colormap(cs)
figure(100)
%yscale("log")
hold on
for oxy=2:16
    jox = squeeze(jox_data(oxy, :));
    sigma_jox = sigma_jox_data(oxy, :);
    c_oxy = squeeze(cOxy_all(2, oxy, :));
    c_oxy_ext = squeeze(cOxy_all(1, oxy, :));
   
    plot(1./c_oxy, 1./jox, "o", "Color",cs(oxy, :))
    % subplot(1, 3, 1)
    % hold on
    % errorbar(r_data(oxy, :), jox", sigma_jox, sigma_jox, "Color",cs(oxy-11, :))
    % subplot(1, 3, 2)
    % hold on
    % plot(r_data(oxy, :), c_oxy, "Color",cs(oxy-11, :))
    % subplot(1, 3, 3)
    % hold on
    % plot(r_data(oxy, :), jox"./c_oxy, "Color",cs(oxy-11, :))
end

%% save results

% save v_max, k_m profiles
save(pp_path + "v_max_profiles_corrected"+string(temp_string)+".mat", "vMax_all")
save(pp_path + "v_max_profiles_corrected_sigma"+string(temp_string)+".mat", "sigma_vMax_all")
save(pp_path + "k_m_profiles_corrected"+string(temp_string)+".mat", "kM_all")
save(pp_path + "k_m_profiles_corrected_sigma"+string(temp_string)+".mat", "sigma_kM_all")

% save corrected oxygen levels
save(pp_path + "cOxy_corr"+string(temp_string)+".mat", "cOxy_all")

% save predicted J_ox gradient
save(pp_path + "jox_pred"+string(temp_string)+".mat", "jox_pred")



%% obtain decay lengths from oxygen gradient
figure(3)
hold on
cs = viridis(16);

R_squared_average = [];
dec_length_collection = [];

fit_start_range = [3]%[1 2 3 4 5 6 7 8 9];
cs = viridis(16);

for fit_start=fit_start_range
    fit_params = [];
    sigma_fit_params = [];
    chi_squareds = [];
    R_squareds = [];
    for oxy=1:numel(precise_oxygen_levels)
        
        norm_log_coxy = log(squeeze(cOxy_all(2, oxy, :))./squeeze(cOxy_all(2, oxy, end)));
        %norm_log_coxy(abs(norm_log_coxy)<1e-8) = 1e-4;
        sigma_norm_log_coxy = (sigma_precise_oxygen_levels./precise_oxygen_levels)
        p0 = [0.06 -1e-4];
        dp = [0.01 0.001];
        
        y_weights = 1./sigma_norm_log_coxy;
        %fit_start = 7;
        fit_end = 10;
    
        pmin = [1e-4 -1e3];
        pmax = [1e4 1e4];
    
        c = [];
    
        % fit
        max_iter = 100;
       
        [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                         lm(@linear_model, p0,...
                         r_data(oxy, :)', norm_log_coxy, y_weights(fit_start:fit_end), dp,...
                         pmin, pmax, c, fit_start, fit_end, max_iter);
        % !!!optional:
        % use sum of squared residuals instead of built-in method of lm to evaluate X2
        x2_jox = sum(((r_data(oxy,fit_start:fit_end) - linear_model(norm_log_coxy(fit_start:fit_end), p, c)).^2)./norm_log_coxy(fit_start:fit_end));
        X2 = x2_jox;

        if oxy > 11
        
        plot(r_data(oxy, :), norm_log_coxy, "o", Color=cs(oxy, :))
        plot(r_data(oxy, :), linear_model(r_data(oxy, :), p, c), Color=cs(oxy, :), LineWidth=1.5)
        
        end
        fit_params = [fit_params, p];
        sigma_fit_params = [sigma_fit_params, sigma_p];
        chi_squareds = [chi_squareds X2];
        R_squareds = [R_squareds R_sq];
    end
    
    decay_lengths = 1./abs(fit_params(1, :));
    sigma_decay_lengths = decay_lengths.*(sigma_fit_params(1, :)./fit_params(1, :));
    R_squared_average = [R_squared_average, mean(R_squareds)];
    size(decay_lengths)
    dec_length_collection = [dec_length_collection; decay_lengths];
end

save(pp_path + "decayLength_coxyEmp" + temp_string + ".mat", "decay_lengths")
save(pp_path + "sigma_decayLength_coxyEmp" + temp_string + ".mat", "sigma_decay_lengths")

%% plot decay length fit parameters for different numbers of included data points

figure(55)
subplot(1, 2, 1)
plot(10-fit_start_range, R_squared_average, "o", "MarkerSize", 15,...
     "LineWidth",1.5)
xlabel("number of data points included in fit")
ylabel("R^2 quality of fit")
title("quality of lin. fit")
% plot(precise_oxygen_levels, fit_params(2, :), "o")
% title("intercept")

cs = viridis(numel(fit_start_range));

subplot(1, 2, 2)
title("decay lengths obtained from fit")
hold on
for fit_start_ind=1:numel(fit_start_range)
    plot(precise_oxygen_levels, dec_length_collection(fit_start_ind, :),...
         "o", "Color",cs(fit_start_ind, :), "MarkerSize",15,...
         "LineWidth",1.5, "DisplayName",string(10-fit_start_range(fit_start_ind))+"fit points")

end
xlabel("c_{oxy} (\mu M)")
ylabel("\lambda_{c_{oxy}} (\mu m)")
legend("Location","southeast")


%% obtain decay lengths from oxygen gradient by fitting sinh(R/lambda)
figure(2)
hold on
cs = viridis(16);

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];
for oxy=1:numel(precise_oxygen_levels)

    % p = [R, A, lambda]
    p0 = [r_data(oxy, end) 1 30];
    dp = [0 0.001 0.001];
    
    y_weights = squeeze(cOxy_all(2, oxy, end))./sigma_precise_oxygen_levels;
    fit_start = 2;
    fit_end = 10;

    pmin = [0 1e-4 1e-4];
    pmax = [1e4 1e4 1e4];

    c = [];

    % fit
    max_iter = 100;
   
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@rd_solver_linear, p0,...
                     r_data(oxy, :)', squeeze(cOxy_all(2, oxy, :))./squeeze(cOxy_all(2, oxy, end)), y_weights(fit_start:fit_end), dp,...
                     pmin, pmax, c, fit_start, fit_end, max_iter);
    % !!!optional:
    % use sum of squared residuals instead of built-in method of lm to evaluate X2
    %x2_jox = sum(((r_data(oxy,fit_start:fit_end) - linear_model(norm_log_coxy(fit_start:fit_end), p, c)).^2)./norm_log_coxy(fit_start:fit_end));
    %X2 = x2_jox;
    
    if oxy<4
    plot(r_data(oxy, :), squeeze(cOxy_all(2, oxy, :))./squeeze(cOxy_all(2, oxy, end)), "o", Color=cs(oxy, :))
    plot(r_data(oxy, :), rd_solver_linear(r_data(oxy, :), p, c), Color=cs(oxy, :), LineWidth=1.5)
    end

    fit_params = [fit_params, p];
    sigma_fit_params = [sigma_fit_params, sigma_p];
    chi_squareds = [chi_squareds X2];
    R_squareds = [R_squareds R_sq];
end

decay_lengths = fit_params(3, :)
sigma_decay_lengths = decay_lengths.*(sigma_fit_params(3, :)./fit_params(3, :));

save(pp_path + "sinhDecayLength_coxy" + temp_string + ".mat", "decay_lengths")
save(pp_path + "sigma_sinhDecayLength_coxy" + temp_string + ".mat", "sigma_decay_lengths")

%%
subplot(1, 2, 1)
plot(precise_oxygen_levels, fit_params(3, :), "o")
title("decay length \lambda")

subplot(1, 2, 2)
plot(precise_oxygen_levels, fit_params(2, :), "o")
title("amplitude A")

%% obtain decay lengths from jox gradient by fitting sinh(R/lambda)
figure(2)
hold on
cs = viridis(16);

fit_params = [];
sigma_fit_params = [];
chi_squareds = [];
R_squareds = [];
for oxy=1:numel(precise_oxygen_levels)

    % p = [R, A, lambda]
    p0 = [r_data(oxy, end) 40 30];
    dp = [0 0.001 0.001];
    
    y_weights = squeeze(cOxy_all(2, oxy, end))./sigma_precise_oxygen_levels;
    fit_start = 2;
    fit_end = 10;

    pmin = [0 1e-4 1e-4];
    pmax = [1e4 1e4 1e4];

    c = [];

    % fit
    max_iter = 100;
   
    [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                     lm(@rd_solver_linear, p0,...
                     r_data(oxy, :)', squeeze(jox_pred(2, oxy, 1:end)), y_weights(fit_start:fit_end), dp,...
                     pmin, pmax, c, fit_start, fit_end, max_iter);
    % !!!optional:
    % use sum of squared residuals instead of built-in method of lm to evaluate X2
    %x2_jox = sum(((r_data(oxy,fit_start:fit_end) - linear_model(norm_log_coxy(fit_start:fit_end), p, c)).^2)./norm_log_coxy(fit_start:fit_end));
    %X2 = x2_jox;
    
    if oxy<4
    plot(r_data(oxy, :), squeeze(jox_pred(2, oxy, 1:end)), "o", Color=cs(oxy, :))
    plot(r_data(oxy, :), rd_solver_linear(r_data(oxy, :), p, c), Color=cs(oxy, :), LineWidth=1.5)
    end

    fit_params = [fit_params, p];
    sigma_fit_params = [sigma_fit_params, sigma_p];
    chi_squareds = [chi_squareds X2];
    R_squareds = [R_squareds R_sq];
end

decay_lengths = fit_params(3, :)
sigma_decay_lengths = decay_lengths.*(sigma_fit_params(3, :)./fit_params(3, :));

save(pp_path + "sinhDecayLength_jox" + temp_string + ".mat", "decay_lengths")
save(pp_path + "sigma_sinhDecayLength_jox" + temp_string + ".mat", "sigma_decay_lengths")

%%
subplot(1, 2, 1)
plot(precise_oxygen_levels, fit_params(3, :), "o")
title("decay length \lambda")

subplot(1, 2, 2)
plot(precise_oxygen_levels, fit_params(2, :), "o")
title("amplitude A")

%% get decay lengths directly from data
figure(3)
hold on
cs = viridis(16);

R_squared_average = [];
dec_length_collection = [];

fit_start_range = [3]%[1 2 3 4 5 6 7 8 9];

for fit_start=fit_start_range

    fit_params = [];
    sigma_fit_params = [];
    chi_squareds = [];
    R_squareds = [];
    cs = viridis(16);
    
    for oxy=1:numel(precise_oxygen_levels)
        
        log_jox = log(jox_data(oxy, :));
        sigma_log_jox = sigma_jox_data(oxy, :)./jox_data(oxy, :);
        p0 = [0.06 -1e-4];
        dp = [0.01 0.001];
        
        y_weights = 1./sigma_log_jox;
        %fit_start = 6;
        fit_end = 10;
    
        pmin = [1e-4 -1e3];
        pmax = [1e4 1e4];
    
        c = [];
    
        % fit
        max_iter = 100;
        size(r_data(oxy, :))
        size(log_jox)
        size(y_weights)
        [p,X2,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
                         lm(@linear_model, p0,...
                         r_data(oxy, :)', log_jox', y_weights(fit_start:fit_end)', dp,...
                         pmin, pmax, c, fit_start, fit_end, max_iter);
        % !!!optional:
        % use sum of squared residuals instead of built-in method of lm to evaluate X2
        % x2_jox = sum(((r_data(oxy,fit_start:fit_end) - linear_model(log_coxy(fit_start:fit_end), p, c)).^2)./log_coxy(fit_start:fit_end));
        % X2 = x2_jox;
        
        plot(r_data(oxy, :), log_jox, "o", Color=cs(oxy, :))
        plot(r_data(oxy, :), linear_model(r_data(oxy, :), p, c), Color=cs(oxy, :), LineWidth=1.5)
    
        fit_params = [fit_params, p];
        sigma_fit_params = [sigma_fit_params, sigma_p];
        chi_squareds = [chi_squareds X2];
        R_squareds = [R_squareds R_sq];
    end
    
    decay_lengths = 1./abs(fit_params(1, :));
    sigma_decay_lengths = decay_lengths.*(sigma_fit_params(1, :)./fit_params(1, :));
    R_squared_average = [R_squared_average, mean(R_squareds)];

    dec_length_collection = [dec_length_collection; decay_lengths];
end

save(pp_path + "decayLength_joxEmp" + temp_string + ".mat", "decay_lengths")
save(pp_path + "sigma_decayLength_joxEmp" + temp_string + ".mat", "sigma_decay_lengths")
%% plot quality of fit for different numbers of included data points

figure(1)
subplot(1, 2, 1)

plot(10-fit_start_range, R_squared_average, "o", "MarkerSize", 15,...
     "LineWidth",1.5)
cs = viridis(numel(fit_start_range));
xlabel("number of data points included in fit")
ylabel("R^2 quality of fit")
title("quality of lin. fit")

subplot(1, 2, 2)
title("decay lengths obtained from fit")
hold on
for fit_start_ind=1:numel(fit_start_range)
    plot(precise_oxygen_levels, dec_length_collection(fit_start_ind, :),...
         "o", "Color",cs(fit_start_ind, :), "MarkerSize",15,...
         "LineWidth",1.5, "DisplayName",string(10-fit_start_range(fit_start_ind))+"fit points")
    if fit_start_ind==3
        plot(precise_oxygen_levels, dec_length_collection(fit_start_ind, :),...
             "LineWidth",1.5, "Color",cs(fit_start_ind, :), "HandleVisibility","off")
    end
end
xlabel("c_{oxy} (\mu M)")
ylabel("\lambda_{c_{oxy}} (\mu m)")
legend("Location","southeast")


