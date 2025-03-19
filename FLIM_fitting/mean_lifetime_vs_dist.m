clear all;

filename_lifetime=dir('*_mean_lifetime_in_rings.mat');

lifetime_all=[];
dist_all=[];


for i=1:length(filename_lifetime)
   
    load(filename_lifetime(i).name);
    
    lifetime_all=[lifetime_all;lifetime_mean_mat_mean];
    dist_all=[dist_all;dist_mito_cell_dist_mat_mean];
    
end

%%

errorbar(nanmean(dist_all),nanmean(lifetime_all),std(lifetime_all)./sqrt(size(lifetime_all,1)),std(lifetime_all)./sqrt(size(lifetime_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Mean NADH lifetime (ns)')
set(gca,'FontSize',15)

ax=gca;
ax.LineWidth=1.5;

