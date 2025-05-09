%% Script for ripping data from multiple 4-plot figs and replotting to compare
% Place script in an analysis folder.
% Results from FLIM_batch_fitting.m should be placed in labeled sub-folders
clear all; close all;
pth = [pwd '\'];
outpth = [pth 'ComparisonPlots\'];
[a,b] = mkdir(outpth);
D = dir(pth); D(1:2)=[];
for i = 1:length(D) dirind(i)=D(i).isdir; end
D = D(dirind);
% Remove 'CompareisonPlots'
for i = 1:length(D)
    if strcmp(D(i).name,'ComparisonPlots') remv = i; end
end
D(remv) = [];

r=[1 0 0]; g=[0 1 0]; b=[0 0 1]; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1]; p=[160 32 240]./255; c = [0 1 1];
col = [b;o;b;o;k;p;c];
markers = {'-.','-.','-','-','x','x'};

%% Rip data from plots
Segs = {'joint','mito','cyto'};
ind = 1:length(D);
% ind = [1 2];

for seg = 2:3
    for i = 1:length(D)
        close all
        uiopen([pth D(i).name '\fits_Pos0_SingleMasks_' Segs{seg} '.fig'],1)
        hs=findall(gcf,'type','axes');
        
        Nirr{i} = hs(2).Children;
        Nirr{i} = [Nirr{i}(end).XData; Nirr{i}(end).YData]'; % 'end' because sometimes misfits are plotted
        Firr{i} = hs(1).Children;
        Firr{i} = [Firr{i}(end).XData; Firr{i}(end).YData]';
        
        
        Nf{i} = hs(6).Children;
        Nf{i} = [Nf{i}(end).XData; Nf{i}(end).YData]';
        Ff{i} = hs(5).Children;
        Ff{i} = [Ff{i}(end).XData; Ff{i}(end).YData]';
        
        labs{i} = D(i).name;
    end
    
    h=figure;
    for i = ind
        plot(Nirr{i}(:,1),Nirr{i}(:,2),markers{i},'color',col(i,:),'linewidth',1.5); hold on;
    end
    axis tight;
    legend(labs(ind),'location','best')
    xlabel('time (h)'); ylabel('Nirr')
    title(Segs{seg})
    saveas(gcf,[outpth Segs{seg} '_Nirr.jpg'])
    saveas(gcf,[outpth Segs{seg} '_Nirr.fig'])
    close(h)
    
    h=figure;
    for i = ind
        plot(Firr{i}(:,1),Firr{i}(:,2),markers{i},'color',col(i,:),'linewidth',1.5); hold on;
    end
    axis tight;
    legend(labs(ind),'location','best')
    xlabel('time (h)'); ylabel('Firr')
    title(Segs{seg})
    saveas(gcf,[outpth Segs{seg} '_Firr.jpg'])
    saveas(gcf,[outpth Segs{seg} '_Firr.fig'])
    close(h)
    
    h=figure;
    for i = ind
        plot(Nf{i}(:,1),Nf{i}(:,2),markers{i},'color',col(i,:),'linewidth',1.5); hold on;
    end
    axis tight;
    legend(labs(ind),'location','best')
    xlabel('time (h)'); ylabel('Nfrac')
    title(Segs{seg})
    saveas(gcf,[outpth Segs{seg} '_Nf.jpg'])
    saveas(gcf,[outpth Segs{seg} '_Nf.fig'])
    close(h)
    
    h=figure;
    for i = ind
        plot(Ff{i}(:,1),Ff{i}(:,2),markers{i},'color',col(i,:),'linewidth',1.5); hold on;
    end
    axis tight;
    legend(labs(ind),'location','best')
    xlabel('time (h)'); ylabel('Ffrac')
    title(Segs{seg})
    saveas(gcf,[outpth Segs{seg} '_Ff.jpg'])
    saveas(gcf,[outpth Segs{seg} '_Ff.fig'])
    close(h)
end
