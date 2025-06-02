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
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4], [4.4, 5.2]};

temp_resolved = false;
data_path = "/home/mx/mitoFluxGradients/data/";
pp_path = '../data/EXP_published/postprocessing_results/';
plot_path = '../data/EXP_published/plots/';
% temp_resolved = true;
% data_path = "/home/mx/mitoFluxGradients/data/";
% pp_path = '../data/EXP_temperatureResolved/postprocessing_results/';
% plot_path = '../data/EXP_temperatureResolved/plots/';

load(data_path + 'exp_solubilities.mat')
load(data_path + 'exp_temperatures.mat')
%%
temp_ind = 1;

if temp_resolved==true
    temp_string = '_T'+string(temperature(temp_ind))+'C'
else
    temp_ind = 4;
    temp_string = '';
end
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
sigma_precise_oxygen_levels = (213.5/20.946).*sigma_precise_oxygen_levels;
sigma_precise_oxygen_levels(sigma_precise_oxygen_levels==0) = mean(sigma_precise_oxygen_levels);
end

%% simple plot of J_ox gradients for different outside oxygen

% subplot with all oxygen levels included
fig = figure('Renderer', 'painters', 'Position', [10 10 600 400]);
cs = flip(viridis(20));
colormap(cs)

hold on;
for k=1:numel(oxygen_ranges)
    
    dist_all=oxygen_ranges_data{k}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{k}.jox_cell_kn_all;

    errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(jox_cell_kn_all(:, 2:end), 'omitnan'),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));

    xlabel('distance to oocyte center (\mu m)');
    ylabel('inferred J_{ox} (\mu M/s)');
    
    %ylim([-20, 80])

    set(gca,'FontSize',15);
end

xlim([4, 36])


cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);

savefig(string(plot_path)+'Jox_vs_dist_emp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_emp'+string(temp_string)+'.png')


% subplot with 3 maximally different levels

fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 400]);
cs = flip(viridis(20));
colormap(cs)

hold on;
for k=[1, 2, numel(oxygen_ranges)]
    
    dist_all=oxygen_ranges_data{k}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{k}.jox_cell_kn_all;

    errorbar(mean(dist_all(:, 2:end), 'omitnan'),mean(jox_cell_kn_all(:, 2:end), 'omitnan'),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             std(jox_cell_kn_all(:, 2:end), 'omitnan')./sqrt(size(jox_cell_kn_all(:, 2:end),1)),...
             'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));

    xlabel('distance to oocyte center (\mu m)');
    ylabel('inferred J_{ox} (\mu M/s)');
    
    %ylim([-20, 80])

    set(gca,'FontSize',15);
end

xlim([4, 36])
%ylim([-25, 140])

cb = colorbar;
clim([1.5 numel(precise_oxygen_levels)+0.5])
title(cb, 'c^* (\mu M)', 'Interpreter', 'tex')
cb.Ticks = linspace(1, numel(precise_oxygen_levels), numel(precise_oxygen_levels));
cb.TickLabels = round(precise_oxygen_levels(1:1:end), 2);

savefig(string(plot_path)+'Jox_vs_dist_emp_maxDiff'+string(temp_string)+'.fig')
saveas(fig2, string(plot_path)+'Jox_vs_dist_emp_maxDiff'+string(temp_string)+'.png')

%% plot Jox as a function of subcellular oxygen for all rings

% choose all or only selection of all distances (rings) to plot
plot_selectedDist = false;
dists = [2, 4, 6, 8, 10];

fig = figure('Renderer', 'painters', 'Position', [10 10 600 400]);
% save corrected oxygen levels
load(string(pp_path)+'cOxy_corr'+string(temp_string)+'.mat')

% subplot(1, 2, 2)
hold on;
cs = flip(plasma(15));
colormap(cs(1:10, 1:end))

xlabel('inferred c(r) (\mu M)');
ylabel('inferred J_{ox}(r, c) (\mu M/s)');
xlim([-5 55])
ylim([-4 135])

set(gca,'FontSize',15);

cb = colorbar;
clim([0.5 10.5])
title(cb, 'r (\mu m)', 'Interpreter', 'tex')
cb.Ticks = linspace(1, 10, 10);
cb.TickLabels = round(mean(dist_all, 1, 'omitnan'), 2);

if plot_selectedDist==false
    for k=2:10
        for j=1:numel(oxygen_ranges)
            dist_all=oxygen_ranges_data{j}.dist_all;
            jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
            errorbar(cOxy_all(2, j, k),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                     std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                     'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));
         
            
        end
    end

    box on
    axes('Position',[.6 .7 .2 .2])
    hold on
    for k=2:10
        for j=1:7
            dist_all=oxygen_ranges_data{j}.dist_all;
            jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
            errorbar(squeeze(cOxy_all(2, j, k)),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                 std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                 'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));
        end    
    end
else
    for k=dists
        for j=1:numel(oxygen_ranges)
            dist_all=oxygen_ranges_data{j}.dist_all;
            jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
            errorbar(cOxy_all(2, j, k),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                     std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                     'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));
        end
    end

    box on
    axes('Position',[.6 .7 .2 .2])
    hold on
    for k=dists
        for j=1:7
            dist_all=oxygen_ranges_data{j}.dist_all;
            jox_cell_kn_all=oxygen_ranges_data{j}.jox_cell_kn_all;
            errorbar(squeeze(cOxy_all(2, j, k)),mean(jox_cell_kn_all(:, k), 'omitnan'),...
                 std(jox_cell_kn_all(:, k), 'omitnan')./sqrt(size(jox_cell_kn_all(:, k),1)),...
                 'o', 'MarkerSize',10,'LineWidth',1.5, 'Color',cs(k, :));
        end    
    end
end



if plot_selectedDist==false
    savefig(string(plot_path)+'Jox_vs_coxy_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'Jox_vs_coxy_emp'+string(temp_string)+'.png')
else
    savefig(string(plot_path)+'Jox_vs_coxy_selected_emp'+string(temp_string)+'.fig')
    saveas(fig, string(plot_path)+'Jox_vs_coxy_selected_emp'+string(temp_string)+'.png')
end
%% calculate mean squared distance (+ ringw-wise covariance, kolmogorow-smirnov test) between oxy levels 

covar = zeros(numel(oxygen_ranges)-1, 10);
kstest2_results = zeros(numel(oxygen_ranges)-1, 10);
kstest2_p = zeros(numel(oxygen_ranges)-1, 10);
klDiv = zeros(numel(oxygen_ranges)-1, 10);
msds = zeros(numel(oxygen_ranges), 0);
corr_matrices = {};
p_vals = {};

for k=1:numel(oxygen_ranges)-1
    k
    current = oxygen_ranges_data{k}.jox_cell_kn_all;
    next = oxygen_ranges_data{k+1}.jox_cell_kn_all;
    trunc_size = min(size(current, 1), size(next, 1));

    msds(k) = mean((mean(current, 'omitnan') - mean(next, 'omitnan')).^2);
    
    [rho,pval] = corr(current(1:trunc_size, :), next(1:trunc_size, :), 'Type', 'Kendall', 'rows', 'complete');
    corr_matrices{k} = rho;
    p_vals{k} = pval;

    for j=1:min(size(current, 2), size(next, 2))
        [kstest2_results(k, j), kstest2_p(k, j)] = kstest2(current(:, j), next(:, j));

        dev_current = std(current(:, j), 'omitnan')./sqrt(size(current(:, j),1));
        dev_next = std(next(:, j), 'omitnan')./sqrt(size(next(:, j),1));

        cov_mat = cov(current(1:trunc_size, j), next(1:trunc_size, j), 'omitrows');
        covar(k, j) = (cov_mat(1, 2) + cov_mat(2, 1))/2;
        covar(k, j) = covar(k, j)/(dev_current*dev_next);
        
    end
end

%% create list of correlation coefficients corresponding to oxygen levels

mean_rho = zeros(numel(oxygen_ranges), 0);
mean_pval = zeros(numel(oxygen_ranges), 0);
for i=1:numel(oxygen_ranges)-1
    mean_rho(i) = mean(diag(corr_matrices{i}), 'all');
    mean_pval(i) = mean(diag(p_vals{i}), 'all');
end

%% plot for comparing gradients at different oxy levels

newcolors = viridis(16);
fig = figure('Renderer', 'painters', 'Position', [10 10 400 300]);
% plot of relative error of J_ox measurement as a function of outside
% oxygen level

rel_err = zeros(numel(precise_oxygen_levels), 10);
for k=1:numel(oxygen_ranges)
    dist_all=oxygen_ranges_data{k}.dist_all;
    jox_cell_kn_all=oxygen_ranges_data{k}.jox_cell_kn_all;
    jox = mean(jox_cell_kn_all, 'omitnan');
    sigma_jox = std(jox_cell_kn_all, 'omitnan')./sqrt(size(jox_cell_kn_all,1));
    rel_err(k, :) = sigma_jox./jox;
end

subplot(2, 2, 1);
hold on;
plot(precise_oxygen_levels, max(rel_err, [], 2), 'o', 'LineWidth',1.5, ...
    'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('max(\sigma_{J_{ox}}/J_{ox})');
set(gca,'FontSize',12);
xlim([0 12])

% plot of correlation coefficient between data set at current
% with data set of next oxygen level
subplot(2, 2, 2);
hold on;
plot(precise_oxygen_levels(1:end-1), mean_rho, 'o', 'LineWidth',1.5,...
    'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('Kendalls tau correlation coefficient');
set(gca,'FontSize',12);
xlim([0 12])

% plot of mean squared distance to next higher J_ox data set
subplot(2, 2, 3);
hold on;
plot(precise_oxygen_levels(1:end-1), msds, 'o', 'LineWidth',1.5,...
    'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('mean squared distance to next set');
set(gca,'FontSize',12);
xlim([0 12])

% plot of KS test results (test, whether current and next J_ox data set
% where sampled from the same prob. distr.
subplot(2, 2, 4);
hold on;
plot(precise_oxygen_levels(1:end-1), mean(kstest2_results, 2),...
    'o', 'LineWidth',1.5, 'DisplayName','rejected h0','MarkerSize',10)

plot(precise_oxygen_levels(1:end-1), mean(kstest2_p, 2), 'o', 'LineWidth',1.5,...
    'Color','red', 'DisplayName','confidence level p', 'MarkerSize',10)
xlabel('oxygen concentration (\mu M)');
ylabel('% of rejected h0 hypothesis (KS test)');
set(gca,'FontSize',12);
xlim([0 12])
legend()

savefig(string(plot_path)+'Jox_vs_dist_empComp'+string(temp_string)+'.fig')
saveas(fig, string(plot_path)+'Jox_vs_dist_empComp'+string(temp_string)+'.png')