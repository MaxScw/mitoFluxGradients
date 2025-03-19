clear all;

filename_cellular_kn=dir('*_jox_in_rings_cellular_kn.mat');

int_all=[];
dist_all=[];
bound_all=[];
long_all=[];
jox_cell_kn_all=[];
short_all=[];


for i=1:length(filename_cellular_kn)
   
    load(filename_cellular_kn(i).name);
    
    int_all=[int_all;irr_mito_cell_dist_mat_mean];
    dist_all=[dist_all;dist_mito_cell_dist_mat_mean];
    bound_all=[bound_all;bound_ratio_mito_nadh];
    long_all=[long_all;long_mito_nadh];
    jox_cell_kn_all=[jox_cell_kn_all;jox];
    short_all=[short_all;short_mito_nadh];
    
end

jox_ring_kn_all=[];
kn_ring_all=[];
long_ring_all=[];

filename_ring_kn=dir('*_jox_in_rings_ring_kn_with_lifetime.mat');


for i=1:length(filename_ring_kn)
   
    load(filename_ring_kn(i).name);
    
    jox_ring_kn_all=[jox_ring_kn_all;jox];
    kn_ring_all=[kn_ring_all;kn_ring];
    long_ring_all=[long_ring_all;long_ring];
    
end

%% FLIM par

subplot(2,2,1)

errorbar(nanmean(dist_all),nanmean(int_all)./mean(int_all(:,end)),std(int_all)./sqrt(size(int_all,1))./mean(int_all(:,end)),std(int_all)./sqrt(size(int_all,1))./mean(int_all(:,end)),'or','MarkerSize',10,'LineWidth',1.5)
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

errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;

subplot(2,2,3)

errorbar(nanmean(dist_all),nanmean(long_all),std(long_all)./sqrt(size(long_all,1)),std(long_all)./sqrt(size(long_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('NADH long lifetime (\tau_{l}) (ns)')
set(gca,'FontSize',15)

%ylim([1.92 2.22])

ax=gca;
ax.LineWidth=2.0;

subplot(2,2,4)

errorbar(nanmean(dist_all),nanmean(short_all),std(short_all)./sqrt(size(short_all,1)),std(short_all)./sqrt(size(short_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('NADH short lifetime (\tau_{s}) (ns)')
set(gca,'FontSize',15)

%ylim([0.4 0.5])

ax=gca;
ax.LineWidth=2.0;

%% Predited Jox

errorbar(nanmean(dist_all),nanmean(jox_cell_kn_all),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Predicted J_{ox} (uM/s)')
set(gca,'FontSize',23)

%ylim([45 85])

ax=gca;
ax.LineWidth=2.5;
%% plot kn

errorbar(nanmean(dist_all),nanmean(kn_ring_all),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
% hold on
% plot(nanmean(dist_all),nanmean(jox_ring_kn_all),'ob')
xlabel('Distance to oocyte center (um)')
ylabel('Equilibrium NADH bound ratio \beta_{eq}')
set(gca,'FontSize',15)

%ylim([0.95 1.2])

ax=gca;
ax.LineWidth=1.5;

%% jox from cell kn vs jox from ring kn

errorbar(nanmean(dist_all),nanmean(jox_cell_kn_all),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
errorbar(nanmean(dist_all),nanmean(jox_ring_kn_all),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),'ob','MarkerSize',10,'LineWidth',1.5)
legend('cellular \beta_{eq}','ring \beta_{eq}')
xlabel('Distance to oocyte center (um)')
ylabel('Predicted J_{ox} (uM/s)')
set(gca,'FontSize',15)

%ylim([45 85])

ax=gca;
ax.LineWidth=1.5;
