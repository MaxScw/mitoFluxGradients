%% read-in of data (not neccessary if read-in data was loaded in cell above)
dpath = 'temp_28C';

for oxy_ind=1:numel(oxygen_ranges)
    oxygen_ranges{oxy_ind}
    filename_cellular_kn=dir([dpath, '/', oxygen_ranges{oxy_ind},'/', '*_jox_in_rings_cellular_kn.mat']);
    filename_irr_decay_dist_in_rings=dir([dpath, '/', oxygen_ranges{oxy_ind},'/', '*_irr_decay_dist_in_rings.mat']);
    oxy.dist_all=[];
    oxy.jox_cell_kn_all=[];
    oxy.o2_levels=[];
    oxy.sigma_o2_levels=[];

    for i=1:length(filename_cellular_kn)
        % filename_cellular_kn(i).name
        load([dpath, '/', oxygen_ranges{oxy_ind},'/', filename_cellular_kn(i).name]);
        
        load([dpath, '/', oxygen_ranges{oxy_ind},'/', filename_irr_decay_dist_in_rings(i).name]);

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

save('plot_data_multiple_oxy_ranges_T28C.mat', 'oxygen_ranges_data');