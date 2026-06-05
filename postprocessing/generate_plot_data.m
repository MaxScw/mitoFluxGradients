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
                   [2.0, 2.8], [2.8, 3.6], [3.6, 4.4],... 
                   [4.4, 5.2]};

load('../data/exp_solubilities.mat')
load('../data/exp_temperatures.mat')

% temp_resolved = false;
% pp_path = '../data/EXP_published/postprocessing_results/';
% plot_path = '../data/EXP_published/plots/';
% dpath = "../data/EXP_published/FLIM_fitting_results/temp_37C";

temp_resolved = true;
pp_path = '/home/mx/mitoFluxGradients/data/EXP_temperatureResolved/postprocessing_results/';
plot_path = '/home/mx/mitoFluxGradients/data/EXP_temperatureResolved/plots';
temp_ind = 5;
if temp_resolved==true
    temp_string = '_T'+string(temperature(temp_ind))+'C';
    dpath = '../data/EXP_temperatureResolved/FLIM_fitting_results/temp_'+string(temperature(temp_ind))+'C';
else 
    temp_string = '';
end



%% read-in of data

for oxy_ind=1:numel(oxygen_ranges)
    oxygen_ranges{oxy_ind};
   
    filename_cellular_kn=dir(dpath + '/' + oxygen_ranges{oxy_ind} + '/' + '*_jox_in_rings_cellular_kn.mat');
    filename_irr_decay_dist_in_rings=dir(dpath + '/' + oxygen_ranges{oxy_ind} + '/' + '*_irr_decay_dist_in_rings.mat');
    oxy.dist_all=[];
    oxy.jox_cell_kn_all=[];
    oxy.o2_levels=[];
    oxy.sigma_o2_levels=[];
    oxy.nadhf=[];
    oxy.nadhb=[];
 

    for i=1:length(filename_cellular_kn)
        % filename_cellular_kn(i).name
        load(dpath + '/' + oxygen_ranges{oxy_ind} + '/' + filename_cellular_kn(i).name);
        
        load(dpath + '/' + oxygen_ranges{oxy_ind} + '/' + filename_irr_decay_dist_in_rings(i).name);

        temp_02_levels = [];
        for fls=1:numel(flim_struct)
            temp_02_levels = [temp_02_levels, flim_struct(fls).oxygen];
        end
        temp_02_levels = temp_02_levels(temp_02_levels>num_oxygen_ranges{oxy_ind}(1));
        temp_02_levels = temp_02_levels(temp_02_levels<num_oxygen_ranges{oxy_ind}(2));
   
        
        oxy.o2_levels = [oxy.o2_levels, mean(temp_02_levels)];
        oxy.sigma_o2_levels = [oxy.sigma_o2_levels, std(temp_02_levels)];
       

        oxy.dist_all=[oxy.dist_all;dist_mito_cell_dist_mat_mean];
        
        oxy.jox_cell_kn_all=[oxy.jox_cell_kn_all;jox];
        oxy.nadhf=[oxy.nadhf;nadhf];
        oxy.nadhb=[oxy.nadhb;nadhf.*(bound_ratio_mito_nadh)];
        %oxy.background=[oxy.background;]

        % if isnan(mean(temp_02_levels))
        %     oxy.dist_all=[oxy.dist_all;NaN];
        %     oxy.jox_cell_kn_all=[oxy.jox_cell_kn_all;NaN];
        % else
        %     oxy.dist_all=[oxy.dist_all;dist_mito_cell_dist_mat_mean];
        %     oxy.jox_cell_kn_all=[oxy.jox_cell_kn_all;jox];
        % end
    end
    
    oxy.o2_levels = mean(oxy.o2_levels, 'omitnan');
    oxy.sigma_o2_levels = mean(oxy.sigma_o2_levels, 'omitnan');
    oxygen_ranges_data{oxy_ind}=oxy;
end

save(string(pp_path)+'plot_data_multiple_oxy_ranges'+string(temp_string)+'.mat', 'oxygen_ranges_data');