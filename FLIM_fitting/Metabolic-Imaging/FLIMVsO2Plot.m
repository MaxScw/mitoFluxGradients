function FLIMVsO2Plot(fitpath,RestBool,O2rng,O2InfoFile)
% Load FLIM decay struct, take O2 concentrations from 'O2FromRuthDecay'
% (saved in exp path of 'fitpath'), and plot FLIM params as a function of
% O2 concentration.
% Inputs:
% -fitpath: path to the fit FILE, not the acquisition folder
% -RestBool: Boolean of whether to plot the FLIM params for the O2 restore.
% -O2rng: optional O2 range to subsample data if desired
% -O2InfoFile: path to the 'O2Info.mat' file. If not specified, assume it
%   is in the same acquisition folder as 'fitpath'

% clear all; close all;
% fitpath = 'C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\2016-11-10 1-cell\s1_a1_O2drop\fits_Pos0_SingleMasks.mat';
% RestBool = 0;
% O2rng = [10^-4 5];

% Versions:
% 2018-06-16: Updated to take mito and cyto segment data.

load(fitpath)
% If not specified, assume that O2Info.mat file was saved to the same
% path as the FLIM acquisition ('fitpath')
if ~exist('O2InfoFile')|O2InfoFile==-1
    O2InfoFile = [UpOneDir(fitpath) '\O2Info.mat'];
end
load(O2InfoFile)

%%
timeunit = 3600; % 60 gives minutes, 3600 for hours, 86400 for days
switch timeunit
    case 1;
        tlab = 'Time (s)';
    case 60;
        tlab = 'Time (min)';
    case 3600;
        tlab = 'Time (h)';
    case 86400;
        tlab = 'Time (days)';
end
% Versions:
% 2016-12-28: changed Drop and Rest ranges to times
% 2016-12-05: Altered to look in root directory for O2info.mat. Also
%    updated plotting with GoodPlotyyTicks.m

%% Group data in one cell. Also construct a cell that can be written in spreadsheet form
AllDecays = [];
labels = {};


tm = [];
chs = [];
% Get t and ch from deay structures instead of relying on nameinds at
% this point
decays_fits_struct(cellfun('isempty',decays_fits_struct)) = [];
RuthExc = [];
for j = 1:size(decays_fits_struct,1)
    slashes = strfind(fitpath,'\');
    [T,ch,Z] = MultiD_Parse_FName(decays_fits_struct{j}.filename);
    if strcmp(ch,'Ruth')
        RuthExc = [RuthExc j];
        continue;
    end
    t = decays_fits_struct{j}.timestp;
    
    tm = [tm; t];
    chs = [chs; {ch}];
    label = [fitpath(slashes(end)+1:end-4) '_' ch '_t' num2str(T) '_z' num2str(Z)];
    labels = [labels; label];
end
currL = length(decays_fits_struct);
decays_fits_struct(RuthExc) = [];
AllDecays = [AllDecays; decays_fits_struct];

% Check if there were multiple segments, 'joint, mito, and cyto' channels, in that order
% If so, create separate plots for each seg with a label. Otherwise, just 
% plot the one segment.
mitobool = 0;
for i = 1:length(AllDecays)
    if ~isempty(find(AllDecays{i}.decay(:,2)))
        mitobool = 1;
    end
end
if mitobool
    segind = 1:size(AllDecays{1}.decay,2);
else
    segind = 1;
end


L = length(AllDecays);
for i = 1:L
    fitres = AllDecays{i}.fit_result;
    AllParams{i} = fitres;
    Chisqred(i,:) = AllDecays{i}.Chi_sq;
    irrs(i,:) = AllDecays{i}.irr;
    numphot(i,:) = sum(AllDecays{i}.decay);
end

% Indices of 4D params matrix:
% ind1 and ind2 are 5x3 fit result matrix
% ind3 is segment (joint, mito, cyto)
% ind4 is decay element number (=sdt number)
for i = 1:length(AllDecays)
    ParamsArr(:,:,:,i) = AllParams{i};
end

for i = 1:L
    %   Plot arrays
    f(i,:) = squeeze(ParamsArr(4,1,:,i))';
    shifts(i,:) = squeeze(ParamsArr(1,1,:,i))';
    tau1(i,:) = squeeze(ParamsArr(3,1,:,i))';
    tau2(i,:) = tau1(i,:).*squeeze(ParamsArr(5,1,:,i))';
end

Nind = strcmp(chs,'NADH');
Find = strcmp(chs,'FAD');


%% Convert from time to O2. Construct O2 array corresponding to tm
for seg = segind
%     seg = 2; % Only plot response of mito channel.
    for i = 1:length(tm)
        dists = abs(tm(i) - O2T);
        ind = find(dists==min(dists)); ind = ind(1);
        O2tm(i) = O2(ind);
    end
    
    % Indices
    % Identify bad fits
    misfits = [];
    misfits = f(:,seg)==0|f(:,seg)==1; %misfits = misfits';
    NDropInd = Nind&~misfits&tm>DropTime&tm<RestTime;
    FDropInd = Find&~misfits&tm>DropTime&tm<RestTime;
    N1Dropind = find(NDropInd); N1Dropind = N1Dropind(1); F1Dropind = find(FDropInd); F1Dropind = F1Dropind(1);
    if RestBool
        NRestInd = Nind&~misfits&tm>RestTime&tm<EndTime;
        FRestInd = Find&~misfits&tm>RestTime&tm<EndTime;
        N1Restind = find(NRestInd); N1Restind = N1Restind(1); F1Restind = find(FRestInd); F1Restind = F1Restind(1);
    else
        NRestInd = zeros(length(Nind),1); FRestInd = zeros(length(Nind),1);
        N1Restind = []; F1Restind = [];
    end
    
    if ~exist('O2rng')|O2rng==-1 O2rng = [min(O2tm(NDropInd|NRestInd|FDropInd|FRestInd)) max(O2tm(NDropInd|NRestInd|FDropInd|FRestInd))]; end
    
    %% If no arguments entered, plot in one big plot.
    % if nargin<2
    fh = gcf;
    clf(fh)
    set(fh,'Name','Parameters Comparison','Units','normalized','Position',[.1 .1 .8 .8],'Color','w'); hold on;
    
    
    
    
    % FBOUNDS
    subplot(2,2,3,'position',[0.1300 0.110 0.3347 0.3712]); hold on;
    title('Population Fractions Bound','fontsize',12)
    grid on;
    % Plot FLIM params for O2 Drop and O2 Restore in different symbols
    %%%%%% ONLY NADH AND FAD CASES CODED. CODE OTHERS IF NEED BE
    if ~isempty(find(Nind))&isempty(find(Find))
        %     [ax,h1,h2] =plotyy(tm(Nind),irrs(Nind),tm(Nind),numphot(Nind));
        %     set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
        %     set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    elseif ~isempty(find(Find))&isempty(find(Nind))
        %     [ax2,h1,h2] =plotyy(tm(Find),irrs(Find),tm(Find),numphot(Find));
        %     set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
        %     set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
        %     xlim([min(tm) max(tm)])
    elseif ~isempty(find(Nind))&~isempty(find(Find))
        [ax,h1,h2] =plotyy(O2tm(NDropInd),f(NDropInd,seg),...
            O2tm(FDropInd),1-f(FDropInd,seg));
        ax(1).XLim = [O2rng]; ax(2).XLim = [O2rng];
        ax(1).YMinorGrid='off'; ax(1).XMinorGrid='off';
        ax(2).YMinorGrid='off'; ax(2).XMinorGrid='off';
        ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
        set(h1,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','b');
        set(h2,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','g');
        hold(ax(1),'on');
        plot(ax(1),O2tm(find(NRestInd)),f(find(NRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','b');
        ax(1).YLim = [min(f(NDropInd|NRestInd,seg))  max(f(NDropInd|NRestInd,seg))];
        hold(ax(2),'on');
        plot(ax(2),O2tm(find(FRestInd)),1-f(find(FRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','g');
        ax(2).YLim = [min(1-f(FDropInd|FRestInd,seg)) max(1-f(FDropInd|FRestInd,seg))];
        ax(1).FontSize = 11; ax(2).FontSize = 11;
        ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
        ax(1).XScale='log'; ax(2).XScale='log';
        % Set reasonable ticks
        GoodPlotyyTicks(ax,O2tm,f(NDropInd|NRestInd,seg),1-f(FDropInd|FRestInd,seg))
        % Mark first point so we can check hysteresis
        plot(ax(1),O2tm(N1Dropind),f(N1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Dropind),1-f(F1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(1),O2tm(N1Restind),f(N1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Restind),1-f(F1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
    end
    set(ax(1),'fontsize',11)
    set(ax(2),'fontsize',11)
    xlabel('O2 %')
    xlhand = get(gca,'xlabel');
    
    % LONG LIFETIMES
    subplot(2,2,1,'position',[0.1300 0.57 0.3347 0.3712]); hold on;
    subplot(2,2,1)
    
    %%%%%% ONLY NADH AND FAD CASES CODED. CODE OTHERS IF NEED BE
    if ~isempty(find(Nind))&isempty(find(Find))
    elseif ~isempty(find(Find))&isempty(find(Nind))
    elseif ~isempty(find(Nind))&~isempty(find(Find))
        [ax,h1,h2] =plotyy(O2tm(NDropInd),tau1(NDropInd,seg),...
            O2tm(FDropInd),tau1(FDropInd,seg));
        ax(1).XLim = [O2rng]; ax(2).XLim = [O2rng];
        ax(1).YMinorGrid='off'; ax(1).XMinorGrid='off';
        ax(2).YMinorGrid='off'; ax(2).XMinorGrid='off';
        ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
        set(h1,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','b');
        set(h2,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','g');
        hold(ax(1),'on');
        plot(ax(1),O2tm(find(NRestInd)),tau1(find(NRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','b');
        ax(1).YLim = [min(tau1(NDropInd|NRestInd,seg))  max(tau1(NDropInd|NRestInd,seg))];
        hold(ax(2),'on');
        plot(ax(2),O2tm(find(FRestInd)),tau1(find(FRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','g');
        ax(2).YLim = [min(tau1(FDropInd|FRestInd,seg)) max(tau1(FDropInd|FRestInd,seg))];
        ax(1).FontSize = 11; ax(2).FontSize = 11;
        ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
        ax(1).XScale='log'; ax(2).XScale='log';
        % Set reasonable ticks
        GoodPlotyyTicks(ax,O2tm,tau1(NDropInd|NRestInd,seg),tau1(FDropInd|FRestInd,seg))
        % Mark first point so we can check hysteresis
        plot(ax(1),O2tm(N1Dropind),tau1(N1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Dropind),tau1(F1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(1),O2tm(N1Restind),tau1(N1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Restind),tau1(F1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
    end
    set(ax(1),'fontsize',11)
    set(ax(2),'fontsize',11)
    xlhand = get(gca,'xlabel');
    grid on;
    title('Long Lifetimes (ns)','fontsize',12)

    
    % SHORT LIFETIMES
    subplot(2,2,2,'position',[0.5703 0.57 0.3347 0.3712]); hold on;
    title('Short Lifetimes (ns)','fontsize',12)
    grid on;
    %%%%%% ONLY NADH AND FAD CASES CODED. CODE OTHERS IF NEED BE
    if ~isempty(find(Nind))&isempty(find(Find))
    elseif ~isempty(find(Find))&isempty(find(Nind))
    elseif ~isempty(find(Nind))&~isempty(find(Find))
        [ax,h1,h2] =plotyy(O2tm(NDropInd),tau2(NDropInd,seg),...
            O2tm(FDropInd),tau2(FDropInd,seg));
        ax(1).XLim = [O2rng]; ax(2).XLim = [O2rng];
        ax(1).YMinorGrid='off'; ax(1).XMinorGrid='off';
        ax(2).YMinorGrid='off'; ax(2).XMinorGrid='off';
        ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
        set(h1,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','b');
        set(h2,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','g');
        hold(ax(1),'on');
        plot(ax(1),O2tm(find(NRestInd)),tau2(find(NRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','b');
        ax(1).YLim = [min(tau2(NDropInd|NRestInd,seg))  max(tau2(NDropInd|NRestInd,seg))];
        hold(ax(2),'on');
        plot(ax(2),O2tm(find(FRestInd)),tau2(find(FRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','g');
        ax(2).YLim = [min(tau2(FDropInd|FRestInd,seg)) max(tau2(FDropInd|FRestInd,seg))];
        ax(1).FontSize = 11; ax(2).FontSize = 11;
        ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
        ax(1).XScale='log'; ax(2).XScale='log';
        % Set reasonable ticks
        GoodPlotyyTicks(ax,O2tm,tau2(NDropInd|NRestInd,seg),tau2(FDropInd|FRestInd,seg))
        % Mark first point so we can check hysteresis
        plot(ax(1),O2tm(N1Dropind),tau2(N1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Dropind),tau2(F1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(1),O2tm(N1Restind),tau2(N1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Restind),tau2(F1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        
    end
    set(ax(1),'fontsize',11)
    set(ax(2),'fontsize',11)
    xlhand = get(gca,'xlabel');
    
    % IRRADIANCE AND PHOTON COUNTS
    subplot(2,2,4,'position',[0.5703 0.110 0.3347 0.3712]); hold on;
    title('Irradiances (photons/area)','fontsize',12)
    grid on;
    %%%%%% ONLY NADH AND FAD CASES CODED. CODE OTHERS IF NEED BE
    if ~isempty(find(Nind))&isempty(find(Find))
    elseif ~isempty(find(Find))&isempty(find(Nind))
    elseif ~isempty(find(Nind))&~isempty(find(Find))
        [ax,h1,h2] =plotyy(O2tm(NDropInd),irrs(NDropInd,seg),...
            O2tm(FDropInd),irrs(FDropInd,seg));
        ax(1).XLim = [O2rng]; ax(2).XLim = [O2rng];
        ax(1).YMinorGrid='off'; ax(1).XMinorGrid='off';
        ax(2).YMinorGrid='off'; ax(2).XMinorGrid='off';
        ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
        set(h1,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','b');
        set(h2,'LineStyle','-','LineWidth',2,'Marker','v','Markersize',11,'Color','g');
        hold(ax(1),'on');
        plot(ax(1),O2tm(find(NRestInd)),irrs(find(NRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','b');
        ax(1).YLim = [min(irrs(NDropInd|NRestInd,seg))  max(irrs(NDropInd|NRestInd,seg))];
        hold(ax(2),'on');
        plot(ax(2),O2tm(find(FRestInd)),irrs(find(FRestInd),seg),'LineStyle','-','LineWidth',2,'Marker','o','Markersize',11,'Color','g');
        ax(2).YLim = [min(irrs(FDropInd|FRestInd,seg)) max(irrs(FDropInd|FRestInd,seg))];
        ax(1).FontSize = 11; ax(2).FontSize = 11;
        ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
        ax(1).XScale='log'; ax(2).XScale='log';
        
        % Set reasonable ticks
        GoodPlotyyTicks(ax,O2tm,irrs(NDropInd|NRestInd,seg),irrs(FDropInd|FRestInd,seg))
        % Mark first point so we can check hysteresis
        plot(ax(1),O2tm(N1Dropind),irrs(N1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Dropind),irrs(F1Dropind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(1),O2tm(N1Restind),irrs(N1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        plot(ax(2),O2tm(F1Restind),irrs(F1Restind,seg),'LineStyle','-','LineWidth',2,'Marker','.','Markersize',11,'Color','r');
        
    end
    % ax(1).YLabel.String =  'NADH Irr';
    % ax(2).YLabel.String =  'FAD Irr';
    grid on;
    set(gca,'fontsize',11)
    xlabel('O2 %')
    % legend(h2,'Photons','location','ne')
    
    % Figure out which case we're dealing with.
    if mitobool
        switch seg
            case 1
                SegLab = '_joint';
            case 2
                SegLab = '_mito';
            case 3
                SegLab = '_cyto';
        end
    else
        SegLab = '';
    end
    set(gcf,'name',[labels{1} SegLab])
    
    % SAVE GRAPHS
    set(gcf,'PaperPositionMode','auto')
    if RestBool
        saveas(gcf,[fitpath(1:end-4) SegLab '_O2Plots_wRest.jpg'],'jpg')
        saveas(gcf,[fitpath(1:end-4) SegLab '_O2Plots_wRest.fig'],'fig')
    else
        saveas(gcf,[fitpath(1:end-4) SegLab '_O2Plots.jpg'],'jpg')
        saveas(gcf,[fitpath(1:end-4) SegLab '_O2Plots.fig'],'fig')
    end
    
end
