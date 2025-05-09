function FLIMBatchParamPlot(names,wrtname,SprBool)
% Plots params from a saved data set from FLIM analysis (FittingGUI_ver3_4)
% names - cell of file paths to fit files, or a string of a single file
% wrtname - optional path to a specific write path
% SprBool - if 1, a spreadsheet is also created with the parameters

% clear all;
% names{1} = 'Z:\Lab\Tim\2017-01-19 Emre collab 3\s1_a1_CLPP\fits_Pos0_mask1_fxshft';
% names{1} = 'C:\Users\Tim\Documents\Academic - Research\Data\fits_Pos0_mask1.mat';
% wrtname = 'C:\Users\Tim\Documents\Academic - Research\Data\2017-05-19 Batch 2\s1_a1\testplot';
% SprBool = 1;

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
% 2018-03-02: Re-write plotting, adding a segment dimension to plot arrays.
% 2018-02-20: Add provision to make separate plots for joint, mito, and,
%    cyto channels
% 2017-07-14: Add 'UserChan' condition that just plots intensity
% 2016-12-05: Put extensive tick creation into a new function: 'GoodPlotyyTicks.m'
% 2016-02-22: Change plotting by a lot. Make double plots to better scale
%    all the data.
% 2015-10-26: % Get t and ch from deay structures instead of relying on
%    nameinds at this point
% 2015-10-10: Revised to stitch trajectories together from multiple decays
% 2015-01-20: Get rid of CalInt and the beads. Obsolete.
% 2014-11-07: Change around plots to display shifts and Chi2 on same plot,
%    then display irradiance in [2,2] subplot.
%   -Also added time(min) for x axis instead of sample labels. Taken from
%    timestamps in nameinds
%   -Also color coding NADH and FAD channels blue and green. Got rid of
%    Anisotropy options because I won't use that for a while.

if ~iscell(names)
    name = names;
    clear names
    names{1} = name;
    clear name
end

%% Group data in one cell. Also construct a cell that can be written in spreadsheet form
AllDecays = [];
labels = {};
% if SprBool
%     % For spreadsheets (one for means, one for std's)
%     spr{1,1} = 'Sample'; spr_std{1,1} = 'Sample';
%     spr{1,2} = 'NADirr'; spr_std{1,2} = 'NADirr';
%     spr{1,3} = 'FADirr'; spr_std{1,3} = 'FADirr';
%     spr{1,4} = 'RedoxRat'; spr_std{1,4} = 'RedoxRat';
%     spr{1,5} = 'NADtau1'; spr_std{1,5} = 'NADtau1';
%     spr{1,6} = 'NADtau2'; spr_std{1,6} = 'NADtau2';
%     spr{1,7} = 'NADfracbound'; spr_std{1,7} = 'NADfracbound';
%     spr{1,8} = 'FADtau1'; spr_std{1,8} = 'FADtau1';
%     spr{1,9} = 'FADtau2'; spr_std{1,9} = 'FADtau2';
%     spr{1,10} = 'FADfracbound'; spr_std{1,10} = 'FADfracbound';
%     spr{1,11} = 'NADtaumean'; spr_std{1,11} = 'NADtaumean';
%     spr{1,12} = 'FADtaumean'; spr_std{1,12} = 'FADtaumean';
% end

tm = [];
chs = [];
UserChBool = 0;
for i=1:length(names)
    load(names{i});
    % Get t and ch from deay structures instead of relying on nameinds at
    % this point.
    decays_fits_struct(cellfun('isempty',decays_fits_struct)) = [];
    % Also remove any frames that are 'Ruth' frames
    RuthExc = [];
    for j = 1:size(decays_fits_struct,1)
        slashes = strfind(names{i},'\');
        [T,ch,Z] = MultiD_Parse_FName(decays_fits_struct{j}.filename);
        if strcmp(ch,'Ruth')
            RuthExc = [RuthExc j];
            continue;
        end
        % Detect whether any user channels are present. If so, don't try to
        % plot any FLIM stuff because it's not there.
        if strcmp(ch,'UserChan') UserChBool = 1; end
        t = decays_fits_struct{j}.timestp;
        
        tm = [tm; t];
        chs = [chs; {ch}];
        %         label = [names{i}(slashes(end)+1:end-4) '_' ch '_t' num2str(T) '_z' num2str(Z)];
        %         labels = [labels; label];
        label = [names{i}(slashes(end)+1:end-4)];
        labels = [labels; label];
    end
    currL = length(decays_fits_struct);
    decays_fits_struct(RuthExc) = [];
    AllDecays = [AllDecays; decays_fits_struct];
    clear decays_fits_struct;
end

tm = (tm - min(tm))*86400/timeunit;
% If there is only one element, and it's t=0, it causes an error in the
% xlim:
if tm==0 tm = [0 0.01]; end

Nind = strcmp(chs,'NADH'); NADHBool = 0;
Find = strcmp(chs,'FAD'); FADBool = 0;
Uind = strcmp(chs,'UserChan');
if ~isempty(find(Nind)) NADHBool = 1; end
if ~isempty(find(Find)) FADBool = 1; end

% if ~UserChBool

% Check if there were multiple segments, 'joint, mito, and cyto' channels, in that order
% If so, create separate plots for each seg with a label. Otherwise, just 
% plot the one segment.
mitobool = 0;
for i = 1:length(AllDecays)
    if size(AllDecays{i}.decay,2)>1 & ~isempty(find(AllDecays{i}.decay(:,2)))
        mitobool = 1;
    end
end
if mitobool
    segind = 1:size(AllDecays{1}.decay,2);
else
    segind = 1;
end


if ~UserChBool
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
        tau2(i,:) = squeeze(ParamsArr(5,1,:,i))';
    end
    
    
    %% Plot in one big plot.
    for seg = segind
        
        % Skip segments that didn't have any content
        if isempty(find(~isnan(irrs(:,seg))))
            continue;
        end
        fh = gcf;
        clf(fh)
        set(fh,'Name','Parameters Comparison','Units','normalized','Position',[.1 .1 .8 .8],'Color','w'); hold on;
        
        % Identify bad fits
        misfits = [];
        misfits = f(:,seg)==0|f(:,seg)==1; %misfits = misfits';
        
        % If no points present, you get plot errors below. Just get the
        % first Nind of Find index and set it equal to 1
        Nind2 = Nind&~misfits;
        Find2 = Find&~misfits;
        if isempty(find(Nind2))&NADHBool
            ind = find(Nind); Nind2(ind(1)) = 1;
        end
        if isempty(find(Find2))&FADBool
            ind = find(Find); Find2(ind(1)) = 1;
        end
        
        
        % SHIFTS and CHI^2's
        subplot(2,2,1,'position',[0.1300 0.57 0.3347 0.3712]); hold on;
        title('Shifts (o) and \chi^2s (.)','fontsize',12)
        grid on;
        if NADHBool&~FADBool
            [ax,h1,h2] =plotyy(tm(Nind2),shifts(Nind2,seg),tm(Nind2),Chisqred(Nind2,seg));
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
        elseif FADBool&~NADHBool
            [ax,h1,h2] =plotyy(tm(Find2),shifts(Find2,seg),tm(Find2),Chisqred(Find2,seg));
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            xlim([min(tm) max(tm)])
        elseif NADHBool&FADBool
            [ax,h1,h2] =plotyy(tm(Nind2),shifts(Nind2,seg),tm(Nind2),Chisqred(Nind2,seg));
            ax(1).XLim = [min(tm) max(tm)];
            ax(2).XLim = [min(tm) max(tm)];
            ax(1).YColor = [0 0 0];
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            hold(ax(1),'on');
            plot(ax(1),tm(Find2),shifts(Find2,seg),'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            ax(1).YLim(1) =  min(shifts(:,seg)); ax(1).YLim(2) =  max(shifts(:,seg));
            hold(ax(2),'on');
            plot(ax(2),tm(Find2),Chisqred(Find2,seg),'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            ylim(ax(2),[min(Chisqred(:,seg))*.95 max(Chisqred(:,seg))*1.05])
            ax(1).FontSize = 11; ax(2).FontSize = 11;
            ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
            % Set reasonable ticks
            ndig = ceil(log10(abs(max(Chisqred(:,seg))-min(Chisqred(:,seg)))));
            ytk1 = round(min(shifts(:,seg)./10^(ndig-1))); ytk1 = ytk1*10^(ndig-1); ytk1 = ytk1(1); ytk2 = round(max(shifts(:,seg)./10^(ndig-1))); ytk2 = ytk2*10^(ndig-1); ytk2 = ytk2(1);
            set(ax(1),'YTick',ytk1:(ytk2-ytk1)/5:ytk2)
            xlim([min(tm) max(tm)])
            % Set reasonable ticks
            ndig = ceil(log10(abs(max(Chisqred(:,seg))-min(Chisqred(:,seg)))));
            ytk1 = round(min(Chisqred(:,seg)./10^(ndig-1))); ytk1 = ytk1*10^(ndig-1); ytk1 = ytk1(1); ytk2 = round(max(Chisqred(:,seg)./10^(ndig-1))); ytk2 = ytk2*10^(ndig-1); ytk2 = ytk2(1);
            set(ax(2),'YTick',ytk1:(ytk2-ytk1)/5:ytk2)
        end
        ax(1).YLabel.String =  'Shifts';
        ax(2).YLabel.String =  '\chi^2s';
        grid on;
        set(gca,'fontsize',12)
        
        
        % FBOUNDS
        subplot(2,2,2)
        subplot(2,2,2,'position',[0.5703 0.57 0.3347 0.3712]);
        title('Population Fractions Bound','fontsize',12)
        grid on; hold on;
        if NADHBool&~FADBool
            plot(tm(Nind2),f(Nind2,seg),'ob','MarkerSize',11,'LineWidth',2);
        elseif FADBool&~NADHBool
            plot(tm(Find2),f(Find2,seg),'og','MarkerSize',11,'LineWidth',2);
        elseif NADHBool&FADBool
            [ax,h1,h2] =plotyy(tm(Nind2),f(Nind2,seg),tm(Find2),1-f(Find2,seg));
            ax(1).XLim = [min(tm) max(tm)];
            ax(2).XLim = [min(tm) max(tm)];
            ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            ax(1).YLabel.String =  'NADH fbound';
            ax(2).YLabel.String =  'FAD fbound';
            ax(1).FontSize = 11; ax(2).FontSize = 11;
            ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
            %     Set reasonable ticks and limits
            GoodPlotyyTicks(ax,f(Nind2,seg),1-f(Find2,seg));
        end
        
        
        % LIFETIMES
        subplot(2,2,3,'position',[0.1300 0.110 0.3347 0.3712]); hold on;
        title('Lifetimes (ns): short (o) and long (\Delta)','fontsize',12)
        grid on;
        twoexpB = 0;
        if NADHBool&~FADBool
            [ax,h1,h2] =plotyy(tm(Nind2),tau1(Nind2,seg),tm(Nind2),tau2(Nind2,seg));
            set(h1,'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
        elseif FADBool&~NADHBool
            [ax2,h1,h2] =plotyy(tm(Find2),tau1(Find2,seg),tm(Find2),tau2(Find2,seg));
            set(h1,'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            xlim([min(tm) max(tm)])
        elseif NADHBool&FADBool
            % Tau1 is short lifetime, tau2 is long
            [ax,h1,h2] =plotyy(tm(Nind2),tau1(Nind2,seg),tm(Nind2),tau2(Nind2,seg));
            ax(1).XLim = [min(tm) max(tm)];
            ax(2).XLim = [min(tm) max(tm)];
            ax(1).YColor = [0 0 0];
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            hold(ax(1),'on');
            plot(ax(1),tm(Find2),tau1(Find2,seg),'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            ax(1).YLim = [min(tau1(~misfits,seg))*.99 max(tau1(~misfits,seg))*1.01];
            hold(ax(2),'on');
            plot(ax(2),tm(Find2),tau2(Find2,seg),'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            ax(2).YLim = [min(tau2(~misfits,seg))*.99 max(tau2(~misfits,seg))*1.01];
            ax(1).FontSize = 11; ax(2).FontSize = 11;
            ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
            %     Set reasonable ticks and limits
            GoodPlotyyTicks(ax,tau1(~misfits,seg),tau2(~misfits,seg));
        end
        ax(1).YLabel.String =  'Short Lifetimes (o)';
        ax(2).YLabel.String =  'Long Lifetimes (\Delta)';
        xlabel(tlab,'fontsize',12)
        
        
        % IRRADIANCES
        subplot(2,2,4,'position',[0.5703 0.110 0.3347 0.3712]); hold on;
        title('Irradiances (phots/pixel/scan)','fontsize',12)
        grid on;
        if NADHBool&~FADBool
            [ax,h1,h2] =plotyy(tm(Nind),irrs(Nind,seg),tm(Nind),numphot(Nind,seg));
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            ax(1).YLabel.String =  'NADH Irr (o)';
            ax(2).YLabel.String =  'NADH Phots';
        elseif FADBool&~NADHBool
            [ax2,h1,h2] =plotyy(tm(Find),irrs(Find,seg),tm(Find),numphot(Find,seg));
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            xlim([min(tm) max(tm)])
            ax(1).YLabel.String =  'FAD Irr (o)';
            ax(2).YLabel.String =  'FAD Phots';
        elseif NADHBool&FADBool
            [ax,h1,h2] =plotyy(tm(Nind),irrs(Nind,seg),tm(Find),irrs(Find,seg));
            ax(1).XLim = [min(tm) max(tm)];
            ax(2).XLim = [min(tm) max(tm)];
            ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
            set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
            set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
            hold(ax(1),'on');
            plot(ax(1),tm(Nind&misfits),irrs(Nind&misfits,seg),'LineStyle','none','Marker','x','Markersize',11,'MarkerEdgeColor','r','LineWidth',2);
            hold(ax(2),'on');
            plot(ax(2),tm(Find&misfits),irrs(Find&misfits,seg),'LineStyle','none','Marker','x','Markersize',11,'MarkerEdgeColor','r','LineWidth',2);
            %             ax(1).YLim = [min(irrs(Nind,seg))*0.99 max(irrs(Nind,seg))*1.01];
            %             ax(2).YLim = [min(irrs(Find,seg))*0.99 max(irrs(Find,seg))*1.01];
            ax(1).FontSize = 11; ax(2).FontSize = 11;
            ax(1).LabelFontSizeMultiplier=1.2; ax(2).LabelFontSizeMultiplier=1.2;
            %     Set reasonable ticks and limits
            GoodPlotyyTicks(ax,irrs(Nind2,seg),irrs(Find2,seg));
            ax(1).YLabel.String =  'NADH Irr';
            ax(2).YLabel.String =  'FAD Irr';
        end
        grid on;
        
        xlabel(tlab,'fontsize',12)
        
        
        % SAVE GRAPHS
        set(gcf,'PaperPositionMode','auto')
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
        if ~isempty(strfind(names{1},'.mat')) names{1} = names{1}(1:end-4); end
        set(gcf,'name',[labels{1} SegLab])
        if size(names,2)==1
            saveas(gcf,[names{1} SegLab '.jpg'],'jpg')
            saveas(gcf,[names{1} SegLab '.fig'],'fig')
        else
            if exist('wrtname')
                if ~isempty(strfind(wrtname,'.mat')) wrtname = wrtname(1:end-4); end
                saveas(gcf,[wrtname SegLab '.jpg'],'jpg')
                saveas(gcf,[wrtname SegLab '.fig'],'fig')
            else
                saveas(gcf,[names{1} SegLab '_JointPlot.jpg'],'jpg')
                saveas(gcf,[names{1} SegLab '_JointPlot.fig'],'fig')
            end
        end
    end
    
else % If UserChan, only look at intensity
    L = length(AllDecays);
    for i = 1:L
        irrs(i) = AllDecays{i}.irrSc;
        numphot(i) = sum(AllDecays{i}.decay);
    end
end



