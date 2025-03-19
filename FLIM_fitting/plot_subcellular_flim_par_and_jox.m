clear all;

filename_cellular_kn=dir(['./', '*_jox_in_rings_cellular_kn.mat']);

int_all=[];
dist_all=[];
bound_all=[];
long_all=[];
jox_cell_kn_all=[];
short_all=[];

for i=1:length(filename_cellular_kn)
    
    load(['./', filename_cellular_kn(i).name]);
    
    int_all=[int_all;irr_mito_cell_dist_mat_mean];
    dist_all=[dist_all;dist_mito_cell_dist_mat_mean];
    bound_all=[bound_all;bound_ratio_mito_nadh];
    long_all=[long_all;long_mito_nadh];
    jox_cell_kn_all=[jox_cell_kn_all;jox];
    short_all=[short_all;short_mito_nadh];
    
end

%% FLIM par

subplot(2,2,1)

errorbar(mean(dist_all, 'omitnan'),mean(int_all, 'omitnan')./mean(int_all(:,end), 'omitnan'),std(int_all, 'omitnan')./sqrt(size(int_all,1))./mean(int_all(:,end), 'omitnan'),std(int_all, 'omitnan')./sqrt(size(int_all,1))./mean(int_all(:,end), 'omitnan'),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Normalized NADH intensity (I)')
set(gca,'FontSize',15)

%ylim([0.95 1.2])

ax=gca;
ax.LineWidth=2.0;

bound_all(isinf(bound_all))=NaN;

subplot(2,2,2)

errorbar(mean(dist_all, 'omitnan'),mean(bound_all, 'omitnan'),std(bound_all, 'omitnan')./sqrt(size(bound_all,1)),std(bound_all, 'omitnan')./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;

subplot(2,2,3)

errorbar(mean(dist_all, 'omitnan'),mean(long_all, 'omitnan'),std(long_all, 'omitnan')./sqrt(size(long_all,1)),std(long_all, 'omitnan')./sqrt(size(long_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('NADH long lifetime (\tau_{l}) (ns)')
set(gca,'FontSize',15)

%ylim([1.92 2.22])

ax=gca;
ax.LineWidth=2.0;

subplot(2,2,4)

errorbar(mean(dist_all, 'omitnan'),mean(short_all, 'omitnan'),std(short_all, 'omitnan')./sqrt(size(short_all,1)),std(short_all, 'omitnan')./sqrt(size(short_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
%plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('NADH short lifetime (\tau_{s}) (ns)')
set(gca,'FontSize',15)

ylim([0.4 0.5])

ax=gca;
ax.LineWidth=2.0;

%% Predited Jox

errorbar(mean(dist_all, 'omitnan'),mean(jox_cell_kn_all, 'omitnan'),std(jox_cell_kn_all, 'omitnan')./sqrt(size(jox_cell_kn_all,1)),std(jox_cell_kn_all, 'omitnan')./sqrt(size(jox_cell_kn_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Predicted J_{ox} (uM/s)')
set(gca,'FontSize',23)

ylim([45 90])

ax=gca;
ax.LineWidth=2.5;
