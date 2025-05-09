% Running 'FLIM_master_list_and_mdpars' creates a master list and multiD
% parameter matrix. Copy this script and those list/MDpar files to the
% anlaysis folder you want, and input the file names below. Then this
% script will make plots and do analysis on the data.

% NOTE: All of these operations treat each sample (embryo/oocytes) as a
% single data point. That means averaging over time points and z-planes,
% but you can specify what Trng and Zrng you want below.

clear all
cd(pwd) % Place a copy of this script in the folder where you want plots to export.
close all
load('colorblind_colormap.mat');
startupTim % initialization to make plots look prettier

%% USER INPUT
OutFile = ['OutName'];
% Load and concatenate lists
listpaths{1} = 'MasterList.xls';
% listpaths{2} = 'MasterList_2.xls';
parpaths{1} = 'MasterList_pars.mat';
% parpaths{2} = 'MasterList_pars2.mat';
[ml,MDpars,txt] = GroupMasterListPars(listpaths,parpaths);
% Convert FAD long fraction to frac 'engaged', inverse
MDpars(4,1,:,:,2,:,:) = 1-MDpars(4,1,:,:,2,:,:);

% Choose segment
segs = [1]; Segs = {'joint','mito','cyto'};

% Choose which analysis to output
FourPlotBool = 1; xyErrBool = 1;
TtestBool = 1;
Plot3Bool = 1;
BswarmBool = 1;

% Specify plot groupings and marker formatting
Grps = {[1],[2]}; % Group sets in cells to plot together, e.g. Grps = {[2 4],[1 3]}
markers = {'o','o','o','o','o','o','o','o'};
r=[1 0 0]; g=[0 1 0]; b=[117 255 255]./255; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1];
col = colorblind;
% col = [b;g]; % Custom colors (comment to use colorblind colors in previous line)
MarkS = 7;

% Specify time range and z-range. Leave empty to average all elements present
Trng = []; if isempty(Trng) Trng = 1:size(MDpars,4); end
Zrng = []; if isempty(Zrng) Zrng = 1:size(MDpars,6); end

%% Auto-setup groups and labels
% Get the group numbers in an array for easy reference. Get Grp inds (row
% indices) for each group.
% If no data set numbers entered, auto-fill with DS1
if isempty(find(~isnan([ml{:,2}]))) ml(:,2) = {1}; end
GrpArr = [ml{:,2}];
for grp = 1:length(Grps)
    Grpinds{grp} = [];
    for j = 1:length(Grps{grp})
        Grpinds{grp} = [Grpinds{grp} find(GrpArr==Grps{grp}(j))];
    end
    % Get group labels from first of grouped sets
    GrpLabs{grp} = ml{Grpinds{grp}(1),3};
    if isnan(GrpLabs{grp}) GrpLabs{grp}=num2str(grp); end
    Ls(grp) = length(Grpinds{grp}); % Array lengths of each group
end
Chlabs = {'NADH', 'FAD', 'User'};
% Detect which channels were present (not just filled with nan's)
ChTstmps = squeeze(nanmean(nanmean(nanmean(nanmean(MDpars(7,1,:,:,:,:,:),3),4),6),7));
chnums = find(~isnan(ChTstmps))';

npar = size(MDpars,1)-1; nexpo = (npar-2)/2;
% Make labels cells in the same order as how they are stored in the
%elements of MDpar. Then use pord to decide what order you want to plot them in.
if nexpo==2
    % Short labels
    ParLabs = {'shift','A(1-bg)','\tau_1','fraction engaged','\tau_2','intensity'};
    UnitsLabs = {'','',' (ns)','',' (ns)',' (au)'};
    pord = [6 4 3 5]; % Default, don't plot shifts and BG
elseif nexpo==3
    ParLabs = {'shift','A(1-bg)','\tau_1','frac2','\tau_2','frac3','\tau_3','intensity'};
    UnitsLabs = {'','',' (ns)','',' (ns)','',' (ns)',' (au)'};
    pord = [8 3 4 5 6 7];
end


% % (Optional) manual legend labels
% LegNum = 2; % Number of legend items to plot
% GrpLabs = {'Young','Old'};

%% AVERAGING:
% AVERAGES CALCULATED WITH WEIGHTING FACTORS = 1/SIGMA^2
% http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
% MDpars dims: [Param#, mean/std(1,2), ListRow, time point, channel, z-position, segment]
% Fixed pars have 0 std err, causing NaNs in means. Just set those to 1.
MDpars(find(MDpars==0))=1; sz = size(MDpars);

ParZaved = nan([sz(1) sz(2) sz(3) sz(4) sz(5) 1 sz(7)]);
if size(MDpars,6)==1
    ParZaved = MDpars;
else
    % Std error - from variance between z values.
    ParZaved(:,1,:,:,:,1,:) = wtnanmean(MDpars(:,1,:,Trng,:,Zrng,:),1./MDpars(:,2,:,Trng,:,Zrng,:),6);
    ParZaved(:,2,:,:,:,1,:) = sqrt(wtnanvar(MDpars(:,1,:,Trng,:,Zrng,:),1./MDpars(:,2,:,Trng,:,Zrng,:),6))/sqrt(size(MDpars,6));
    %     ParZstdnw = sqrt(wtnanvar(MDpars(:,1,:,Trng,:,Zrng,:),ones(size(MDpars(:,1,:,Trng,:,Zrng,:))),6))/sqrt(size(MDpars,6)); % non-weighted
end

% For final average, consider some time frames might be missing, or we may
% be comparing single frames to each other (e.g. before/after a
% perturbation). So for each row, find out how many time points are ~nans.
ParAved = nan([sz(1) sz(2) sz(3) 1 sz(5) 1 sz(7)]);
[r,c]=find(~isnan(squeeze(ParZaved(7,1,:,:,1,1,1))));
for rw = 1:size(ParZaved,3)
    Trng2 = c(r==rw);
    Trng2 = intersect(Trng,Trng2);
    if length(Trng2)==1
        ParAved(:,:,rw,1,:,:,:) = ParZaved(:,:,rw,Trng2,:,:,:);
    else
        ParAved(:,1,rw,1,:,:,:) = wtnanmean(ParZaved(:,1,rw,Trng2,:,1,:),1./ParZaved(:,2,rw,Trng2,:,1,:),4);
        ParAved(:,2,rw,1,:,:,:) = sqrt(wtnanvar(ParZaved(:,1,rw,Trng2,:,1,:),1./ParZaved(:,2,rw,Trng2,:,1,:),4))/sqrt(length(Trng2));
    end
end


%% Plot averages by group in 'Four-plots'
if FourPlotBool
    for seg = segs
        if nexpo==3 error('Only available for 2-exp'); end
        pord = [ 1 2 6 4 3 5];
        close all;
        figure('position',[200 100 1200 750]) % Widescreen
        %     figure('position',[200 100 1000 1000]) % equal-axis ratios
        FontS = 16;
        for p = 3:npar
            %         for ch = 1:2
            dat = nan(max(Ls),length(Grpinds),size(ParAved,5));
            % Fill in dat colums with group values
            subplot(2,2,p-2)
            for grp = 1:length(Grpinds)
                x = squeeze(ParAved(pord(p),1,Grpinds{grp},1,1,1,seg)); % NADH means
                y = squeeze(ParAved(pord(p),1,Grpinds{grp},1,2,1,seg)); % FAD means
                dx = squeeze(ParAved(pord(p),2,Grpinds{grp},1,1,1,seg)); % NADH std errs
                dy = squeeze(ParAved(pord(p),2,Grpinds{grp},1,2,1,seg)); % FAD std errs
                h = plot(x,y,markers{grp},'markersize',MarkS,'markerfacecolor',col(grp,:),'color',col(grp,:));hold on
                if xyErrBool % Error bars
                    Mx = wtnanmean(x,1./dx); % Weighted distribution mean for group
                    Vx = wtnanvar(x,1./dx); % Weighted std dev
                    Ex = sqrt(Vx./length(x)).*1.96; % 95% conf intervals
                    My = wtnanmean(y,1./dy); % Weighted distribution mean for group
                    Vy = wtnanvar(y,1./dy); % Weighted std dev
                    Ey = sqrt(Vy./length(y)).*1.96; % 95% conf intervals
                    errorbar(Mx,My,Ey,Ey,Ex,Ex,'color',col(grp,:)*.7,'linewidth',2)
                end
            end
            xlabel([Chlabs{1} ' ' ParLabs{pord(p)} UnitsLabs{pord(p)}],'fontsize',FontS*1.1)
            ylabel([Chlabs{2} ' ' ParLabs{pord(p)} UnitsLabs{pord(p)}],'fontsize',FontS*1.1)
            set(gca,'fontsize',FontS); set(gcf,'color','w')
            TightWSpace(gca,.07)
        end
        ResizeSubplots(gcf,[.06 .06 .05 .04],[.08 .08])
        saveas(gcf,[OutFile '_FourPlot_' Segs{seg} '.png']);
    end
end


%% T-tests
if TtestBool
    for seg = segs
        % Define a spreadsheet to easily display all the T-test results
        spr = {'Sample','T-test pass','p value','Mean1','Mean2','Var/mn^2_1','Var/mn^2_2','Relmean','Relbar'}
        % For multiple segments, use 3D table
        % spr(:,:,2) = spr(:,:,1);
        
        % if TtestBool
        for p = 1:npar
            for ch = 1:2
                if ch==1 offs = 0; end
                if ch==2 offs = 8; end
                dat = nan(max(Ls),2,length(Grpinds));
                % Fill in dat colums with group values
                for grp = 1:length(Grpinds)
                    dat(1:length(Grpinds{grp}),1,grp) = squeeze(ParAved(p,1,Grpinds{grp},1,ch,1,seg));
                    dat(1:length(Grpinds{grp}),2,grp) = squeeze(ParAved(p,2,Grpinds{grp},1,ch,1,seg));
                end
                x = dat(:,1,1); y = dat(:,1,2); dx = dat(:,2,1); dy = dat(:,2,2);
                L = length(find(~isnan(x)));
                [h,pval,ci,stats]=ttest2(x,y);
                
                % Mean relative difference (mean of ratios)
                F = y./x; dF = abs(F).*sqrt( (dx./x).^2 + (dy./y).^2);
                
                spr{p+1+offs,2} = h;
                spr{p+1+offs,3} = pval;
                spr{p+1+offs,4} = wtnanmean(x,1./dx); spr{p+1+offs,5} = wtnanmean(y,1./dy);
                spr{p+1+offs,6} = wtnanvar(x,1./dx)/wtnanmean(x,1./dx)^2; spr{p+1+offs,7} = wtnanvar(y,1./dy)/wtnanmean(y,1./dy)^2;
                spr{p+1+offs,8} = wtnanmean(F,1./dF)-1; spr{p+1+offs,9} = 1.96*sqrt(wtnanvar(F,1./dF))/sqrt(L);
                Lab = ParLabs{p}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
                spr{p+1+offs,1} = [Chlabs{ch} '_' Lab];
            end
        end
        delete([OutFile '_Ttests_' Segs{seg} '.xls'])
        xlswrite([OutFile '_Ttests_' Segs{seg} '.xls'],spr)
    end
end


%% 3D Plot
if Plot3Bool
    for seg = segs
        close all
        % Manually pick which dimensions to plot. We could look at T-test values
        % above to automatically pick the 3 best separating parameters (but I
        % haven't done that).
        PlotDims = [ 5 4 5];
        ChDims   = [ 2 1 1];
        meas = []; species = [];
        FontS = 16;
        for grp = 1:length(Grpinds)
            x = squeeze(ParAved(PlotDims(1),1,Grpinds{grp},1,ChDims(1),1,seg));
            y = squeeze(ParAved(PlotDims(2),1,Grpinds{grp},1,ChDims(2),1,seg));
            z = squeeze(ParAved(PlotDims(3),1,Grpinds{grp},1,ChDims(3),1,seg));
            h = plot3(x,y,z,markers{grp},'markersize',8,'markerfacecolor',col(grp,:),'color',col(grp,:));hold on
            
            % Build vars for SVM
            meas = [meas; [x y z]];
            species = [species; ones(length(x),1).*grp];
            
            % Get group names
            GroupNames{grp} = [ml{Grpinds{grp}(1),3} ' (n=' num2str(length(x)) ')'];
        end
        xlabel([Chlabs{ChDims(1)} ' ' ParLabs{PlotDims(1)} UnitsLabs{PlotDims(1)}],'fontsize',FontS*1.1)
        ylabel([Chlabs{ChDims(2)} ' ' ParLabs{PlotDims(2)} UnitsLabs{PlotDims(2)}],'fontsize',FontS*1.1,'position',[1.831 0.358 1.802])
        zlabel([Chlabs{ChDims(3)} ' ' ParLabs{PlotDims(3)} UnitsLabs{PlotDims(3)}],'fontsize',FontS*1.1)
        
        set(gca,'fontsize',FontS,'position',[0.13948      0.15039      0.77372      0.77461])
        set(gcf,'position',[571   183   800   527])
        grid on
        axis tight
        
        % DRAW A PLANE
        % Use SVM to find the optimal plane for
        [SVMModel,ph] = SVM_3dPlot_plane_TS(meas,species,[0,1;1,0],100);
        
        % Manual view and manipulations
        % view([-139.1         29.2])
        % % view([0 0])
        %     rotate(ph,[1 1 -1],1)
        %     rotate(ph,[0 0 1],-8)
        %     rotate(ph,[-1 -1 0],-1);
        
        % rotate(ph,[0 0 1],45)
        % % rotate(ph,[1 0 0],15)
        % ylim([1.87 2.23])
        % zlim([.18 .4])
        
        legend(GroupNames,'location','best')%,'position',[0.2009      0.77114      0.27125      0.15655])
        saveas(gcf,[OutFile '_3DPlot_' Segs{seg} '.png']);
    end
end

%% Plot averages by group in bwarms
if BswarmBool
    
    % Set figure sizes
    x1 = 1; x2 = 2; MidW = .15; CapS = 20; Crop = 0.7; pos = [560   528   413   330]; Spread = 13;
    xs = [1 2];
    FontS = 20;
    
    figure;
    for p = 1:npar
        for ch = 1:2
            dat = nan(max(Ls),length(Grpinds));
            % Fill in dat colums with group values
            for grp = 1:length(Grpinds)
                dat(1:length(Grpinds{grp}),grp) = squeeze(ParAved(p,1,Grpinds{grp},1,ch,1,seg));
            end
            close all; h = UnivarScatter(dat,'Label',GrpLabs,'MarkerFaceColor',col,'MarkerEdgeColor',col,'Whiskers','none'); hold on;
            for grp = 1:length(Grpinds)
                E = ErrorBars95Conf(dat(:,grp)); M = nanmean(dat(:,grp)); S = nanstd(dat(:,grp));
                he=errorbar(xs(grp),M,E,'linewidth',1.5,'color','k');
                he.CapSize = CapS;
                plot([xs(grp)-MidW xs(grp)+MidW],[M M],'-k','linewidth',1.5);
            end
            axis tight;
            xlim([min(xs)-Crop max(xs)+Crop]);
            try
                rnge = max(dat(:))-min([dat(:)]);
                ylim([min([dat(:)])-rnge/20 max(dat(:))+rnge/20]);
            catch end
            set(gca,'fontsize',FontS); set(gcf,'color','w')
            ylabel([Chlabs{ch} ' ' ParLabs{p} UnitsLabs{p}],'fontsize',FontS*1.1)
            set(gca,'PlotBoxAspectRatioMode','auto')
            set(gcf,'position',pos)
            axpos = get(gca,'position'); set(gca,'position',axpos+[0 0 0.06 0])
            my_xticklabels(gca,[1 2],GrpLabs,'fontangle','italic','fontsize',FontS);
            Lab = ParLabs{p}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
            Lab(Lab==' ')=[]; % Remove spaces
            saveas(gcf,[OutFile '_Bswrm_' Chlabs{ch} '_' Lab '_' Segs{seg} '.png']);
        end
    end
end
