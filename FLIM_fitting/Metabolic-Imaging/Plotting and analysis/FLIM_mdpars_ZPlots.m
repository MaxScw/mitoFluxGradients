% Running 'FLIM_master_list_and_mdpars' creates a master list and multiD
% parameter matrix. Copy this script and those list/MDpar files to the
% anlaysis folder you want, and input the file names below. Then this
% script will make z-frame plots.

% Merge B2 and B3 data sets, color NADH blue, FAD green.
% Also do t-test stats to show how often 1st division looks different from
% compaction, and compaction from late blast.

clear all
cd(pwd) % Place a copy of this script in the folder where you want plots to export.
close all
load('colorblind_colormap.mat');
startupTim % initialization to make plots look prettier

%% USER INPUT
OutFile = ['ZPlots_'];
% Load and concatenate lists
listpaths{1} = 'MasterList.xls';
% listpaths{2} = 'MasterList2.xls';
parpaths{1} = 'MasterList_pars.mat'; %'MasterList_pars.mat';
% parpaths{2} = 'MasterList_pars2.mat';
[ml,MDpars,txt] = GroupMasterListPars(listpaths,parpaths);
% Convert FAD long fraction to frac 'engaged', inverse
MDpars(4,1,:,:,2,:,:) = 1-MDpars(4,1,:,:,2,:,:);

% Choose segments to plot
segs = [1]; Segs = {'joint','mito','cyto'};

% Choose ValTypes to plot
% -> ValType: 1 = Average absolute values, 2 = Average absolute change,
%        3 = Average percentage change
ValInds = [1];

% Choose which analysis to output
SubPlotsBool = 1; % All params in one fig w/ subplots
AvePlotsBool = 1; % If you want to plot the average curve in bold, enter '1'
BGShfBool = 0; % Whether to display BG and shift on subplots
IndivPlotsBool = 0; % Individual time plots with larger axes areas.

% Choose whether to interpolate

% Specify plot groupings and marker formatting
Grps = {[1],[2]}; % Group sets in cells to plot together, e.g. Grps = {[2 4],[1 3]}
markers = {'o','o','o','o','o','o','o','o'};
r=[1 0 0]; g=[0 1 0]; b=[0 0 1]; bl=[117 255 255]./255; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1];
col = colorblind; %
% col = [b;g]; % Custom colors (comment to use colorblind colors in previous line)
LineW = 2; % Line width
MarkS = 7; % Marker size
LineTrans = 0.4; % Line transparency

% Set time range.
% 1) [] Default  - plot full time range, including far time points that only
%               have few or only one trajectory.
%               NOTE: average curves at these extremes show large
%               fluctuations due to small numbers of trajectories.
%               Typically, you want to crop at some width where you have
%               plenty of trajectories and the average is well-behaved.
% 2) [Tmin,Tmax] (array of doubles) - set custom time range
TDispRan = [];

% Time frame offset. Enter -1 if you don't want to shift time courses.
% 6 is the 1st division, 7 is 2nd division, 8 is blastocoel
TmFrOffsetColm = -1;

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

% Select appropriate time units. If more that 90 min, use hours
AllTs = MDpars(7,1,:,:,:,:,:,:);
AllTs = MDpars(7,1,[Grpinds{:}],:,:,:,:,:);

ExpT = (max(AllTs(:))-min(AllTs(:)))*86400/60;
% 60 gives minutes, 3600 for hours, 86400 for days
if ExpT < 90
    tlab = 'Time (min)';
    timeunit = 60;
else
    tlab = 'Time (h)';
    timeunit = 3600;
end

ValTypeLabs = {'abs. value','abs. change','% change'};
if BGShfBool ParLabs = {'shift','background','intensity','fraction engaged','\tau_1','\tau_2'}; end

% % (Optional) manual legend labels
% LegNum = 2; % Number of legend items to plot
% GrpLabs = {'Young','Old'};

%% AVERAGING: mainly to average across Z to get single value for each time point

% Don't need this for simple z-plots

%% Time-sync'ing

% Don't need this for simple z-plots

%% Calculate average time curves

% Don't need this for simple z-plots

%% Load different calculates (plot types) into array to loop over

% Don't need this for simple z-plots

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Plot all 8 params on single plot with Subplots

% Don't need this for simple z-plots

%% Z-PLOTS

% Determine if there are multiple z-planes. If so, plot as a function of
% frame. Also create an additional x-axis that shows the physical z
% position (taken from 'AcqSettings.txt')

% MDpars dims: [1=Param#, 2=mean/std(1,2), 3=ListRow, 4=time point, 5=channel, 6=z-position, 7=segment]
% mdpars: [1=Param#, 2=mean/stderr(1,2), 3=time point, 4=position, 5=channel, 6=z-position, 7=embryo number, 8=segment]

Xlbs = {'Frame #','Z (\mum)'};
Xflbs = {'XaxFr','XaxZum'};
% -> Could adapt to pot as function as microns, but this would take a
% little time, so only do if necessary.

for xlb = 1%:2 % 1 plots with 'frames' x-axis. Set to 2 to plot vs microns
    
    vt = 1; % Used value types for ParsCalc, but again, only code that if necessary
    
    % Set figure sizes
    MidW = .15; CapS = 20; Crop = 0.7; pos = [51 113 1300 684];
    Spread = 13; FontS = 16;
    if ~BGShfBool pord = [6 4 3 5]; pos = [51 113 1600 684]; end
    close all;
    
    for seg = 1:size(MDpars,7)
        subind = 1;
        for p = 1:length(pord)
            for ch = 1:2
                subplot(length(chnums),length(pord),p+length(pord)*(ch-1))
                
                for grp = 1:length(Grpinds)
                    % Fill in dat colums with group values
                    dat = permute(MDpars(pord(p),1,Grpinds{grp},:,ch,:,seg,vt),[6 4 3 1 2 5 7]);
                    % [6 3 7 1 2 4 5 8]
                    dat = reshape(dat,size(MDpars,4)*size(MDpars,6),length(Grpinds{grp}));
                    % Frames go down rows, columns are separate masks, if present
                    daterr = permute(MDpars(pord(p),2,Grpinds{grp},:,ch,:,seg,vt),[6 4 3 1 2 5 7]);
                    daterr = reshape(daterr,size(MDpars,4)*size(MDpars,6),length(Grpinds{grp}));
                    
                    if xlb==1
                        Xarr = repmat((1:size(dat,1))',1,length(Grpinds{grp})); % Frames
                    else
                        Xarr = repmat(Zs',size(MDpars,3),length(Grpinds{grp})); % z-plane in microns
                    end
                    % Dims of ParsCalc: [par, ListRow, time(synth), channel, seg]
                    if vt==3 daterr = daterr*100; dat = dat*100; end % Convert to percent
                    if ~isempty(dat(~isnan(dat)))
                        plts = plot(Xarr,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW); hold on;
                        pltsE = errorbar(gca,Xarr,dat,daterr,'color',[col(grp,:) LineTrans],'linestyle','none'); hold on;
                        %                         for pl = 1:length(plts)
                        %                             plts(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                        %                             plts(pl).Color = [col(pl,:) LineTrans];
                        %                             pltsE(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                        %                             pltsE(pl).Color = [col(pl,:) LineTrans];
                        %                             PltLbs{pl} = ['m' num2str(pl)];
                        %                         end
                        
                        for pl = 1:size(dat,2)
                            x = Xarr(:,pl); y = dat(:,pl); yerr = daterr(:,pl);
                            rm = isnan(x); x(rm)=[]; y(rm)=[]; yerr(rm)=[];
                            if ~isempty(x)
                                plts(pl) = plot(x,y,'color',[col(1,:) LineTrans],'linewidth',LineW); hold on;
                                pltsE(pl) = errorbar(x,y,yerr,'color',[col(1,:) LineTrans],'linestyle','none'); hold on;
                                plts(pl).DisplayName = [Chlabs{ch} '_rw' Grpinds{grp}(pl)];
                                plts(pl).Color = [col(grp,:) LineTrans];
                                pltsE(pl).DisplayName = [Chlabs{ch} '_rw' Grpinds{grp}(pl)];
                                pltsE(pl).Color = [col(grp,:) LineTrans];
                                %                                 PltLbs{pl} = ['m' num2str(Grpinds{grp}(pl))];
                            end
                        end
                        
                        h(grp) = plts(1); % For legend
                        
                        %% Don't need average plots YET. But can code if necessary
                        %                         GrpAveX = nanmean(Xarr,2); GrpAveY = squeeze(GrpAves(pord(p),1,:,posnum,ch,1,1,seg,vt));
                        %                         GrpAveErr = squeeze(GrpAves(pord(p),2,:,posnum,ch,1,1,seg,vt));
                        %                         rm = isnan(GrpAveX); GrpAveX(rm)=[]; GrpAveY(rm)=[]; GrpAveErr(rm)=[];
                        %                         %             if AvePlotsBool
                        %                         if vt==3 GrpAveErr = GrpAveErr*100; GrpAveY = GrpAveY*100; end % Convert to percent
                        %                         errorbar(gca,GrpAveX,GrpAveY,GrpAveErr,'color','k','LineWidth',1.5);
                        %                 h(1) = plot(GrpAveX,GrpAveY,'color',col(1,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_Grp' num2str(Grps{1})]);
                        %             end
                    end
                    %                     TightWSpace(gca,.05)
                    %                         xlim(TDispRan)
                end
                set(gca,'fontsize',FontS); set(gcf,'color','w')
                xlabel(Xlbs{xlb});
                
                if vt==1 YValLab = UnitsLabs{pord(p)}; else YValLab = [' ' ValTypeLabs{vt}]; end
                ylab = [Chlabs{ch} ' ' ParLabs{pord(p)} YValLab];
                if pord(p)==4&vt>1 ylab = {[Chlabs{ch} ' ' ParLabs{pord(p)}]; YValLab}; end
                if vt==2 ylab = ['\Delta' Chlabs{ch} ' ' ParLabs{pord(p)}]; end
                ylabel(ylab,'fontsize',FontS);
            end
        end
        set(gcf,'position',pos)
        ResizeSubplots(gcf,[.04 .05 .02 .04],[.05 .07])
        legend(h,GrpLabs,'location','best','fontsize',FontS);
        xlabel(Xlbs{xlb},'position',[.32 -.142],'fontsize',FontS)
        legend(h,GrpLabs,'position',[0.92521    0.0042981     0.073125       0.2288]);
        flab = ValTypeLabs{vt}; flab(strfind(flab,'.'))=[]; flab(strfind(flab,' '))=[];
        saveas(gcf,[OutFile Chlabs{ch} '_' flab '_' Segs{seg} '.png']);
    end
end

