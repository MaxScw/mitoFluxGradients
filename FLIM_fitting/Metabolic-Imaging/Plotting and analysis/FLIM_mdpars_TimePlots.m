% Running 'FLIM_master_list_and_mdpars' creates a master list and multiD
% parameter matrix. Copy this script and those list/MDpar files to the
% anlaysis folder you want, and input the file names below. Then this
% script will make time course plots.

% Merge B2 and B3 data sets, color NADH blue, FAD green.
% Also do t-test stats to show how often 1st division looks different from
% compaction, and compaction from late blast.

clear all
cd(pwd) % Place a copy of this script in the folder where you want plots to export.
close all
load('colorblind_colormap.mat');
startupTim % initialization to make plots look prettier

%% USER INPUT
OutFile = ['TPlots_'];
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
RedoxPlotBool = 1; % Individual time plots with larger axes areas.

% Choose whether to interpolate 

% Specify plot groupings and marker formatting
Grps = {[1]}; % Group sets in cells to plot together, e.g. Grps = {[2 4],[1 3]}
markers = {'o','o','o','o','o','o','o','o'};
r=[1 0 0]; g=[0 1 0]; b=[0 0 1]; bl=[117 255 255]./255; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1];
col = colorblind; %
% col = [b;g]; % Custom colors (comment to use colorblind colors in previous line)
LineW = 1; % Line width
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
% AVERAGES CALCULATED WITH WEIGHTING FACTORS = 1/SIGMA^2
% http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance

% MDpars dims: [1=Param#, 2=mean/std(1,2), 3=ListRow, 4=time point, 5=channel, 6=z-position, 7=segment]

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
% Same dims as MDpars, but 1 element in z.

% % In case we want to look at Z dependence (we don't, generally)
% ParsT = squeeze(ParZaved(:,1,:,:,:,:,:)); % Take par vals only
% % Dims of ParsT: [par, ListRow, time, channel, seg].
%
% % In case we want to look at Z dependence
% ParTaved = wtnanmean(MDpars(:,1,:,Trng,:,Zrng,:),1./MDpars(:,2,:,Trng,:,Zrng,:),4);
% ParTaved(:,2,:,:,:,:,:) = sqrt(1./nansum(1./MDpars(:,2,:,Trng,:,Zrng,:).^2,4)); % Propagated error
% ParsZ = squeeze(ParTaved);

%% Time-sync'ing
% Perform optional time synchronization. Also convert to hours and start
% at 't=0' (not absolute time)
ParsTsnc = ParZaved;
for rw = 1:size(ParsTsnc,3)
    for ch = 1:size(ParsTsnc,5)
        for sg = 1:size(ParsTsnc,7)
            if TmFrOffsetColm>-1
                TmFrOffset = ml{rw,TmFrOffsetColm};
                ParsTsnc(7,1,rw,:,ch,1,sg) = (ParsTsnc(7,1,rw,:,ch,1,sg)-ParsTsnc(7,1,rw,TmFrOffset,ch,1,sg))*86400/timeunit;
                % Note, start time is 1st channel. FAD measurements come
                % slightly after NADH.
            else % if no sync specified, just start each traj at 0
                ParsTsnc(7,1,rw,:,ch,1,sg) = (ParsTsnc(7,1,rw,:,ch,1,sg)-ParsTsnc(7,1,rw,1,ch,1,sg))*86400/timeunit;
            end
        end
    end
end

%% Calculate average time curves
% Consider each time point to be an average of all the data available at
% that time point. Since these embryos are sync'd by division times, an
% example is the '2-cell division' time point would be an average of all
% the embryos where the first division was captured. Use linear
% interpolation to have a consistent, oversampled time axis that is
% shared by all curves.
% Note: trajectories are not extrapolated beyond their original time range.
% Note: if a sample has an missing time point, interp between the adjacent
% points.
% Define a synthetic time axis that will be common to all interpolated
% trajectories
% Dims of ParsTsnc: [par, ListRow, time, channel, seg].
AllTs = ParsTsnc(7,1,[Grpinds{:}],:,:,:,:);
T = min(AllTs(:)):(max(AllTs(:))-min(AllTs(:)))/1000:max(AllTs(:));
if ~exist('TDispRan')|isempty(TDispRan) TDispRan = [min(T) max(T)]; end
T(T<TDispRan(1)|T>TDispRan(2))=[];
sz = size(ParsTsnc);
% Interp params
ParsItrp = nan(sz(1),sz(2),sz(3),size(T,2),sz(5),sz(6),sz(7));
for p = 1:npar
    for rw = 1:size(ParsTsnc,3)
        for ch = 1:size(ParsTsnc,5)
            for sg = 1:size(ParsTsnc,7)
                % Interp params
                x = squeeze(ParsTsnc(7,1,rw,:,ch,1,sg));
                y = squeeze(ParsTsnc(p,1,rw,:,ch,1,sg));
                ystd = squeeze(ParsTsnc(p,2,rw,:,ch,1,sg));
                nanind = isnan(x) | isnan(y);
                x(nanind)=[]; y(nanind)=[]; ystd(nanind)=[];
                % Fill in dat colums with group values
                if ~isempty(x)&length(x)>1 % sometimes a segment is missing, e.g. no mito pixels for dim embryo
                    ParsItrp(7,1,rw,:,ch,1,sg) = T;
                    ParsItrp(p,1,rw,:,ch,1,sg) = interp1(x,y,T);
                    % Also carry weights over to interpolation at discrete
                    % time points for later plotting error bars.
                    ParsItrp(p,2,rw,:,ch,1,sg) = interp1(x,ystd,T);
                end
                %                                 % Test
                %                                 plot(squeeze(ParsItrp(7,rw,:,ch,1,sg)),squeeze(ParsItrp(p,rw,:,ch,1,sg)));hold on;
                %                                 plot(x,y,'.r');
                %                                 1;
            end
        end
    end
end

%% Load different calculates (plot types) into array to loop over
% Have multiple types of plots in this analysis:
% -> ValType: 1 = Average absolute values, 2 = Average absolute change,
%        3 = Average percentage change
ParsCalc = nan([size(ParsItrp) 3]);
ParsCalc(:,:,:,:,:,:,:,1) = ParsItrp; % Abs

% CALCULATIONS
t0 = 1; % Frame of time point to divide time course by.

% Abs change
ParsCalc(:,:,:,:,:,:,:,2) = ParsItrp;
X = ParsItrp(1:6,1,:,:,:,:,:); dX = ParsItrp(1:6,2,:,:,:,:,:);
Y = ParsItrp(1:6,1,:,t0,:,:,:); dY = ParsItrp(1:6,2,:,t0,:,:,:);
F = X-Y; % Function - Relative change
ParsCalc(1:6,1,:,:,:,:,:,2) = F;
ParsCalc(1:6,2,:,:,:,:,:,2) = sqrt( dX.^2 + dY.^2 ); % Prop uncert

% Relative change
ParsCalc(:,:,:,:,:,:,:,3) = ParsItrp;
X = ParsItrp(1:6,1,:,:,:,:,:); dX = ParsItrp(1:6,2,:,:,:,:,:);
Y = ParsItrp(1:6,1,:,t0,:,:,:); dY = ParsItrp(1:6,2,:,t0,:,:,:);
F = X./Y; % Function - Relative change
ParsCalc(1:6,1,:,:,:,:,:,3) = F-1; 
ParsCalc(1:6,2,:,:,:,:,:,3) = abs(F).*sqrt( (dX./X).^2 + (dY./Y).^2); % Prop uncert

% Redox ratio(s)
ParsCalc(:,:,:,:,:,:,:,4) = ParsItrp;
X = ParsItrp(1:6,1,:,:,1,:,:); dX = ParsItrp(1:6,2,:,:,1,:,:);
Y = ParsItrp(1:6,1,:,:,2,:,:); dY = ParsItrp(1:6,2,:,t0,2,:,:);
F = X./Y; % Function - NADH pars/FAD pars
ParsCalc(1:6,1,:,:,1,:,:,4) = F;
ParsCalc(1:6,2,:,:,1,:,:,4) = abs(F).*sqrt( (dX./X).^2 + (dY./Y).^2); % Prop uncert

sz = size(ParsItrp);
GrpAves = nan([sz(1) sz(2) length(Grpinds) sz(4) sz(5) 1 sz(7) 4]);
for grp = 1:length(Grpinds)
    GrpAves(1:6,1,grp,:,:,:,:,:) = wtnanmean(ParsCalc(1:6,1,Grpinds{grp},:,:,:,:,:),1./ParsCalc(1:6,2,Grpinds{grp},:,:,:,:,:),3);
    % GrpAves(1,:,:,:,:) = nanmean(ParsCalc(:,1,:,:,:,:),3); % Optional test against straight mean calc
    GrpAves(1:6,2,grp,:,:,:,:,:) = sqrt(wtnanvar(ParsCalc(1:6,1,Grpinds{grp},:,:,:,:,:),1./ParsCalc(1:6,2,Grpinds{grp},:,:,:,:,:),3))./sqrt(length(Grpinds{grp})); % Std errs
    GrpAves(7,:,grp,:,:,:,:,:,:) = ParsCalc(7,:,Grpinds{grp}(1),:,:,:,:,:); % Times shouldn't get averaged, just transferred
end

% % Test
% plot(squeeze(GrpAves(7,1,:,ch,1,sg)),squeeze(GrpAves(4,1,:,ch,1,sg)));hold on;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Plot all 8 params on single plot with Subplots
if SubPlotsBool
    for vt = ValInds % value type -> different calculations derived in ParsCalc
        % Set figure sizes
        MidW = .15; CapS = 20; Crop = 0.7; pos = [51 113 1300 684];
        Spread = 13; FontS = 16;
        if BGShfBool pord = [1 2 6 3 4 5]; pos = [51 113 1600 684]; end
        close all;
        
        % Define a time gap to plot spaced out, interpolated error bars
        % Find non-nan time indices. Squeeze into a matrix of nanBools,
        % then sum down columns... it works.
        Nntimeinds = find(sum(~isnan(ParsTsnc(7,1,[Grpinds{:}],:,1,1,1)),3));
        Tgap = round(1000/(length(Nntimeinds)-1));
        % add factor change spacing if you want to make them look nicely spaced
        
        for seg = segs
            for p = 1:length(pord)
                for ch = chnums
                    subplot(length(chnums),length(pord),p+length(pord)*(ch-1))
                    for grp = 1:length(Grpinds)
                        % Fill in dat colums with group values
                        % Dims of ParsCalc: [par, ListRow, time(synth), channel, seg]
                        dat = squeeze(ParsCalc(pord(p),1,Grpinds{grp},:,ch,1,seg,vt))';
                        daterr = squeeze(ParsCalc(pord(p),2,Grpinds{grp},:,ch,1,seg,vt));
                        T = squeeze(ParsCalc(7,1,Grpinds{grp},:,ch,1,seg,vt))';
                        
                        if vt==3 daterr = daterr*100; dat = dat*100; end % Convert to percent
                        plts = plot(T,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW); hold on;
                        
                        h(1) = plts(1);
                        GrpAveX = squeeze(GrpAves(7,1,grp,:,ch,1,seg,vt)); GrpAveY = squeeze(GrpAves(pord(p),1,grp,:,ch,1,seg,vt));
                        GrpAveErr = squeeze(GrpAves(pord(p),2,grp,:,ch,1,seg,vt));
                        rm = isnan(GrpAveX); GrpAveX(rm)=[]; GrpAveY(rm)=[]; GrpAveErr(rm)=[];
                        if AvePlotsBool
                            if vt==3 GrpAveErr = GrpAveErr*100; GrpAveY = GrpAveY*100; end % Convert to percent
                            errorbar(gca,GrpAveX(1:Tgap:end),GrpAveY(1:Tgap:end),GrpAveErr(1:Tgap:end),'color',col(grp,:),'LineWidth',1.5);
                        end
                    end
                    
                    TightWSpace(gca,.05)
                    xlim(TDispRan)
                    set(gca,'fontsize',FontS); set(gcf,'color','w')
                    xlabel(tlab);
                    if vt==1 YValLab = UnitsLabs{pord(p)}; else YValLab = [' ' ValTypeLabs{vt}]; end
                    ylab = [Chlabs{ch} ' ' ParLabs{pord(p)} YValLab];
                    if pord(p)==4&vt>1 ylab = {[Chlabs{ch} ' ' ParLabs{pord(p)}]; YValLab}; end
                    if vt==2 ylab = ['\Delta' Chlabs{ch} ' ' ParLabs{pord(p)}]; end
                    ylabel(ylab,'fontsize',FontS);
                end
            end
            set(gcf,'position',pos)
            ResizeSubplots(gcf,[.04 .05 .02 .04],[.05 .07])
            %                 legend(plts,PltLbs,'location','best','fontsize',8);
            xlabel(tlab,'position',[.32 -.142],'fontsize',FontS)
            %                 legend(plts,PltLbs,'position',[0.92521    0.0042981     0.073125       0.2288]);
            flab = ValTypeLabs{vt}; flab(strfind(flab,'.'))=[]; flab(strfind(flab,' '))=[];
            %             if BGsclBool flab = [flab '_BGscl']; end
            %             if BGShfBool flab = [flab '_BGShft']; end
            saveas(gcf,[OutFile '_Tsubplots_' flab '_' Segs{seg} '.png']);
            %             saveas(gcf,[OutFile '_Tsubplots_' flab '_' Segs{seg} '.fig']);
        end
    end
end


%% Individual Time plots
% Kind of always easier to look at all params together in one
% subplots... but here you go
if IndivPlotsBool
    % Set figure sizes
    MidW = .15; CapS = 20; Crop = 0.7; pos = [420   278   560   480];
    Spread = 13; FontS = 16;
    close all;
    
    for vt = ValInds % value type -> different calculations derived in ParsCalc
        %         Define a time gap to plot spaced out, interpolated error bars
        Nntimeinds = find(sum(~isnan(ParsTsnc(7,1,[Grpinds{:}],:,1,1,1)),3));
        Tgap = round(1000/(length(Nntimeinds)-1));
        
        for seg = segs
            for p = 1:length(pord)
                for ch = chnums
                    hf = figure('position',pos);
                    for grp = 1:length(Grpinds)
                        % Fill in dat colums with group values
                        % Dims of ParsCalc: [par, ListRow, time(synth), channel, seg]
                        dat = squeeze(ParsCalc(pord(p),1,Grpinds{grp},:,ch,1,seg,vt))';
                        daterr = squeeze(ParsCalc(pord(p),2,Grpinds{grp},:,ch,1,seg,vt));
                        T = squeeze(ParsCalc(7,1,Grpinds{grp},:,ch,1,seg,vt))';
                        
                        if vt==3 daterr = daterr*100; dat = dat*100; end % Convert to percent
                        plts = plot(T,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW); hold on;
                        
                        h(1) = plts(1);
                        GrpAveX = squeeze(GrpAves(7,1,grp,:,ch,1,seg,vt)); GrpAveY = squeeze(GrpAves(pord(p),1,grp,:,ch,1,seg,vt));
                        GrpAveErr = squeeze(GrpAves(pord(p),2,grp,:,ch,1,seg,vt));
                        rm = isnan(GrpAveX); GrpAveX(rm)=[]; GrpAveY(rm)=[]; GrpAveErr(rm)=[];
                        if AvePlotsBool
                            if vt==3 GrpAveErr = GrpAveErr*100; GrpAveY = GrpAveY*100; end % Convert to percent
                            errorbar(gca,GrpAveX(1:Tgap:end),GrpAveY(1:Tgap:end),GrpAveErr(1:Tgap:end),'color',col(grp,:),'LineWidth',1.5);
                        end
                    end
                    
                    TightWSpace(gca,.05)
                    xlim(TDispRan)
                    set(gca,'fontsize',FontS); set(gcf,'color','w')
                    xlabel(tlab);
                    ylab = [Chlabs{ch} ' ' ParLabs{pord(p)} ' ' ValTypeLabs{vt}];
                    if p==4
                        ylab = {[Chlabs{ch} ' ' ParLabs{pord(p)}]; ValTypeLabs{vt}};
                    end
                    if vt==2 ylab = ['\Delta' Chlabs{ch} ' ' ParLabs{pord(p)}]; end
                    ylabel(ylab,'fontsize',FontS);
                    Lab = ParLabs{pord(p)}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
                    Lab(Lab==' ')=[]; % Remove spaces
                    flab = ValTypeLabs{vt}; flab(strfind(flab,'.'))=[]; flab(strfind(flab,' '))=[];
                    saveas(gcf,[OutFile Chlabs{ch} '_' Lab '_' flab '_' Segs{seg} '.png']);
                    close(hf)
                end
            end
        end
    end
end


%% Redox Time plots

if RedoxPlotBool
    % Set figure sizes
    MidW = .15; CapS = 20; Crop = 0.7; pos = [420   278   560   480];
    Spread = 13; FontS = 16;
    close all;
    
    for vt = 4 % value type -> different calculations derived in ParsCalc
        %         Define a time gap to plot spaced out, interpolated error bars
        Nntimeinds = find(sum(~isnan(ParsTsnc(7,1,[Grpinds{:}],:,1,1,1)),3));
        Tgap = round(1000/(length(Nntimeinds)-1));
        
        ch = 1;
        for seg = segs
            hf = figure('position',pos);
            for grp = 1:length(Grpinds)
                % Fill in dat colums with group values
                % Dims of ParsCalc: [par, ListRow, time(synth), channel, seg]
                dat = squeeze(ParsCalc(6,1,Grpinds{grp},:,ch,1,seg,vt))';
                daterr = squeeze(ParsCalc(6,2,Grpinds{grp},:,ch,1,seg,vt));
                T = squeeze(ParsCalc(7,1,Grpinds{grp},:,ch,1,seg,vt))';
                
                if vt==3 daterr = daterr*100; dat = dat*100; end % Convert to percent
                plts = plot(T,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW); hold on;
                
                h(1) = plts(1);
                GrpAveX = squeeze(GrpAves(7,1,grp,:,ch,1,seg,vt)); GrpAveY = squeeze(GrpAves(6,1,grp,:,ch,1,seg,vt));
                GrpAveErr = squeeze(GrpAves(6,2,grp,:,ch,1,seg,vt));
                rm = isnan(GrpAveX); GrpAveX(rm)=[]; GrpAveY(rm)=[]; GrpAveErr(rm)=[];
                if AvePlotsBool
                    if vt==3 GrpAveErr = GrpAveErr*100; GrpAveY = GrpAveY*100; end % Convert to percent
                    errorbar(gca,GrpAveX(1:Tgap:end),GrpAveY(1:Tgap:end),GrpAveErr(1:Tgap:end),'color',col(grp,:),'LineWidth',1.5);
                end
            end
            
            TightWSpace(gca,.05)
            xlim(TDispRan)
            set(gca,'fontsize',FontS); set(gcf,'color','w')
            xlabel(tlab);
            ylab = ['Redox Ratio'];
            ylabel(ylab,'fontsize',FontS);
            saveas(gcf,[OutFile '_RedoxRatio_' Segs{seg} '.png']);
            close(hf)
        end
    end
end