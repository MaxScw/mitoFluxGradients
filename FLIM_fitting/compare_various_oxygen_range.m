clear all;

load('oxy_0_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(jox_cell_kn_all),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
hold on
% load('oxy_0p5_to_1.mat')
% errorbar(nanmean(dist_all),nanmean(jox_cell_kn_all),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),'ob','MarkerSize',10,'LineWidth',1.5)
% hold on
load('oxy_2_to_5p2.mat')
errorbar(nanmean(dist_all),nanmean(jox_cell_kn_all),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),std(jox_cell_kn_all)./sqrt(size(jox_cell_kn_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
legend('0.26 to 5uM oxygen','20 to 52uM oxygen')
xlabel('Distance to oocyte center (um)')
ylabel('Predicted J_{ox} (uM/s)')
set(gca,'FontSize',23)

%ylim([45 85])

ax=gca;
ax.LineWidth=2.5;
%%

load('oxy_0_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(jox_ring_kn_all),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
hold on
% load('oxy_0p5_to_1.mat')
% errorbar(nanmean(dist_all),nanmean(jox_ring_kn_all),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),'ob','MarkerSize',10,'LineWidth',1.5)
% hold on
load('oxy_2_to_5p2.mat')
errorbar(nanmean(dist_all),nanmean(jox_ring_kn_all),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),std(jox_ring_kn_all)./sqrt(size(jox_ring_kn_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
legend('0.26 to 5uM oxygen','20 to 52uM oxygen')
xlabel('Distance to oocyte center (um)')
ylabel('Predicted J_{ox} (uM/s)')
set(gca,'FontSize',23)

%ylim([45 85])

ax=gca;
ax.LineWidth=2.5;


%%
load('oxy_0p1_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
load('oxy_0_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'om','MarkerSize',10,'LineWidth',1.5)
hold on
errorbar(nanmean(dist_all),nanmean(kn_ring_all),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
legend('1 to 5uM oxygen','0.26 to 5uM oxygen','0.26uM oxygen')
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;


%%
load('oxy_0p1_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(long_all),std(long_all)./sqrt(size(long_all,1)),std(long_all)./sqrt(size(long_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
load('oxy_0_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(long_all),std(long_all)./sqrt(size(long_all,1)),std(long_all)./sqrt(size(long_all,1)),'om','MarkerSize',10,'LineWidth',1.5)
hold on
load('oxy_0_to_0p5_with_long.mat')
errorbar(nanmean(dist_all),nanmean(long_ring_all),std(long_ring_all)./sqrt(size(long_ring_all,1)),std(long_ring_all)./sqrt(size(long_ring_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)

legend('1 to 5uM oxygen','0.26 to 5uM oxygen','0.26uM oxygen')
xlabel('Distance to oocyte center (um)')
ylabel('NADH long lifetime')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;

%%
load('oxy_0_to_0p3.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
errorbar(nanmean(dist_all),nanmean(kn_ring_all),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;

%%
load('oxy_0_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
errorbar(nanmean(dist_all),nanmean(kn_ring_all),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;

%%
load('oxy_0_to_0p3.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
errorbar(nanmean(dist_all),nanmean(kn_ring_all),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;

%%
load('oxy_0_to_0p2.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
errorbar(nanmean(dist_all),nanmean(kn_ring_all),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),std(kn_ring_all)./sqrt(size(kn_ring_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;


%%

load('oxy_0_to_0p5.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'or','MarkerSize',10,'LineWidth',1.5)
hold on
load('oxy_0_to_0p3.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'ob','MarkerSize',10,'LineWidth',1.5)
hold on
load('oxy_0_to_0p2.mat')
errorbar(nanmean(dist_all),nanmean(bound_all),std(bound_all)./sqrt(size(bound_all,1)),std(bound_all)./sqrt(size(bound_all,1)),'ok','MarkerSize',10,'LineWidth',1.5)
xlabel('Distance to oocyte center (um)')
ylabel('NADH bound ratio (\beta)')
set(gca,'FontSize',15)

%ylim([0.55 0.78])

ax=gca;
ax.LineWidth=2.0;
