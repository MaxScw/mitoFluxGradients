clear all;


filename_cellular_kn=dir('.\mito_density_oocytes\*_mitotracker_*int_fitted_in_rings.mat')
filename_mito_num=dir('.\mito_density_oocytes\*_mitotracker_*mito_num_in_rings.mat')

int_all_mitotracker=[];
dist_all_mitotracker=[];
mito_num_all_mitotracker=[];
embr_num_all_mitotracker=[];
ratio_num_all_mitotracker=[];
density_num_all_mitotracker=[];


for i=1:length(filename_cellular_kn)
   
    load(filename_cellular_kn(i).name);
    load(filename_mito_num(i).name);
    
    int_all_mitotracker=[int_all_mitotracker;irr_mito_cell_dist_mat_mean];
    dist_all_mitotracker=[dist_all_mitotracker;dist_mito_cell_dist_mat_mean];
    mito_num_all_mitotracker=[mito_num_all_mitotracker;mito_pixel_num_mean];
    embr_num_all_mitotracker=[embr_num_all_mitotracker;embr_pixel_num_mean];
    ratio_num_all_mitotracker=[ratio_num_all_mitotracker;mito_pixel_num_mean./embr_pixel_num_mean];
    density_num_all_mitotracker=[density_num_all_mitotracker;irr_mito_cell_dist_mat_mean.*mito_pixel_num_mean./embr_pixel_num_mean];

    
end



%%

errorbar(nanmean(dist_all_mitotracker),nanmean(int_all_mitotracker)./mean(int_all_mitotracker(:,end)),std(int_all_mitotracker)./sqrt(size(int_all_mitotracker,1))./mean(int_all_mitotracker(:,end)),std(int_all_mitotracker)./sqrt(size(int_all_mitotracker,1))./mean(int_all_mitotracker(:,end)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Normalized mitotracker intensity')
set(gca,'FontSize',15)

ylim([0.8 2.2])

ax=gca;
ax.LineWidth=1.5;

%%

errorbar(nanmean(dist_all_mitotracker),nanmean(mito_num_all_mitotracker),std(mito_num_all_mitotracker)./sqrt(size(mito_num_all_mitotracker,1)),std(mito_num_all_mitotracker)./sqrt(size(mito_num_all_mitotracker,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Mito pixel number')
set(gca,'FontSize',15)

%ylim([0.8 2.2])

ax=gca;
ax.LineWidth=1.5;

%%

errorbar(nanmean(dist_all_mitotracker),nanmean(embr_num_all_mitotracker),std(embr_num_all_mitotracker)./sqrt(size(embr_num_all_mitotracker,1)),std(embr_num_all_mitotracker)./sqrt(size(embr_num_all_mitotracker,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Oocyte pixel number')
set(gca,'FontSize',15)

%ylim([0.8 2.2])

ax=gca;
ax.LineWidth=1.5;

ax=gca;
ax.LineWidth=1.5;



%%

errorbar(nanmean(dist_all_mitotracker),nanmean(ratio_num_all_mitotracker),std(ratio_num_all_mitotracker)./sqrt(size(ratio_num_all_mitotracker,1)),std(ratio_num_all_mitotracker)./sqrt(size(ratio_num_all_mitotracker,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Mito pixel density')
set(gca,'FontSize',15)

%ylim([0.8 2.2])

ax=gca;
ax.LineWidth=1.5;

ax=gca;
ax.LineWidth=1.5;


%%

errorbar(nanmean(dist_all_mitotracker),nanmean(density_num_all_mitotracker),std(density_num_all_mitotracker)./sqrt(size(density_num_all_mitotracker,1)),std(density_num_all_mitotracker)./sqrt(size(density_num_all_mitotracker,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Mito density')
set(gca,'FontSize',15)

%ylim([0.8 2.2])

ax=gca;
ax.LineWidth=1.5;

ax=gca;
ax.LineWidth=1.5;


