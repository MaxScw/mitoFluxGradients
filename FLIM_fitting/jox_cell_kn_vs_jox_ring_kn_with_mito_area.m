clear all;

filename_cellular_kn=dir('*_jox_in_rings_cellular_kn.mat');
filename_mito_area=dir('*_mito_area_in_rings.mat');

int_all=[];
dist_all=[];
bound_all=[];
long_all=[];
jox_cell_kn_all=[];
short_all=[];
mitoarea_all=[];


for i=1:length(filename_cellular_kn)
   
    load(filename_cellular_kn(i).name);
    load(filename_mito_area(i).name);
    
    int_all=[int_all;irr_mito_cell_dist_mat_mean];
    dist_all=[dist_all;dist_mito_cell_dist_mat_mean];
    bound_all=[bound_all;bound_ratio_mito_nadh];
    long_all=[long_all;long_mito_nadh];
    jox_cell_kn_all=[jox_cell_kn_all;jox];
    short_all=[short_all;short_mito_nadh];
    mitoarea_all=[mitoarea_all;mito_size_mat_mean];
    
end

jox_ring_kn_all=[];
kn_ring_all=[];

filename_ring_kn=dir('*_jox_in_rings_ring_kn.mat');


for i=1:length(filename_ring_kn)
   
    load(filename_ring_kn(i).name);
    
    jox_ring_kn_all=[jox_ring_kn_all;jox];
    kn_ring_all=[kn_ring_all;kn_ring];
    
end

%% integrated flux

flux_all_celluar_kn=mitoarea_all.*jox_cell_kn_all;
flux_all_celluar_kn_mean=nanmean(flux_all_celluar_kn);
integrated_flux_celluar_kn_mean=sum(flux_all_celluar_kn_mean);

flux_all_ring_kn=mitoarea_all.*jox_ring_kn_all;
flux_all_ring_kn_mean=nanmean(flux_all_ring_kn);
integrated_flux_ring_kn_mean=sum(flux_all_ring_kn_mean);

%% averaged flux

flux_all_celluar_kn=mitoarea_all.*jox_cell_kn_all;
sum_flux_all_celluar_kn=sum(flux_all_celluar_kn,2);
mito_area_all_sum=sum(mitoarea_all,2);
averaged_flux_celluar_kn_mean=nanmean(sum_flux_all_celluar_kn./mito_area_all_sum);


flux_all_ring_kn=mitoarea_all.*jox_ring_kn_all;
sum_flux_all_ring_kn=sum(flux_all_ring_kn,2);
mito_area_all_sum=sum(mitoarea_all,2);
averaged_flux_ring_kn_mean=nanmean(sum_flux_all_ring_kn./mito_area_all_sum);

