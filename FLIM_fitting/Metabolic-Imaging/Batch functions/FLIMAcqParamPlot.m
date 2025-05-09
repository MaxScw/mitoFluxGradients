function FLIMAcqParamPlot(acqpath,Excs,ValInds,BGsclBool,BGShfBool)
% Plots params from a saved data set from mdpars (multi-d parameter matrix)
% Adapted from FLIM_mdpars_TimePlots to act on individual positions of
% acquisition. Quicker way to visualize acquisitions that doing plots for
% individual masks.
% Z-scans: Code automatically detects whether acquisition was a time lapse
% or zscan. Makes appropriate plots

% INPUTS:
% -acqpath: duh
% -Excs: exclusions. nx2 array of pos-mask pairs. Can rerun this after
%   initial analysis to exclude bad masks.
% -ValInds: optional index for extra kinds of plotting. Enter array for
%   desired plots
%    -> ValTypes: 1 = Average absolute values, 2 = Average absolute change,
%        3 = Average percentage change 
% -BGsclBool: We've found that the background intensity drifts quite a lot.
%   Calculate background intensity over time with 'Make_Masks_from_IlastikProb.m',
%   then this can load that data to scale (divide) by the background
%   intensity to adjust for intensity drift.

% clear all;
% acqpath = 'C:\Dropbox\data\s1_a1_zscan\';
% % Excs = [0 1];
% ValInds = [1];
% BGsclBool = 0;

if acqpath(end)~='\' acqpath = [acqpath '\']; end
if ~exist('Excs')|Excs==-1 Excs = []; end
if ~exist('ValInds')|ValInds==-1 ValInds = 1; end
if ~exist('BGsclBool')|BGsclBool==-1 BGsclBool = 0; end
if ~exist('BGShfBool')|BGShfBool==-1 BGShfBool = 0; end
% Whether to display BG and shift on subplots
startupTim

load('colorblind_colormap.mat'); % In Metab code folders
colorblind = [colorblind; rand(10,3)]; % Add rands in case of more than 12 masks in one position (rare)
colorblind(5,:)=[]; % Delete black, that's used for averages
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end
% For all multiD_pars files present
MDir = dir([acqpath 'multiD_pars*.mat']);
for mdnum = 1:length(MDir)
% Load mdpars
load([acqpath MDir(mdnum).name]);
MDlab = MDir(mdnum).name(13:end-4);
% Select appropriate time units. If more that 90 min, use hours
AllTs = mdpars(7,1,:,:,:,:,:,:);
ExpT = (max(AllTs(:))-min(AllTs(:)))*86400/60;

% 60 gives minutes, 3600 for hours, 86400 for days
if ExpT < 90
    tlab = 'Time (min)';
    timeunit = 60;
else
    tlab = 'Time (h)';
    timeunit = 3600;
end

% Determine which channels are present
uchans = unique(nameinds(:,4));
if ~isempty(find(strcmp(uchans,'NADH'))) Chind(1) = 1; end
if ~isempty(find(strcmp(uchans,'FAD'))) Chind(2) = 2; end
if ~isempty(find(strcmp(uchans,'UserChan'))) Chind(3) = 3; end
Chind(Chind==0)=[];
ChLabs = {'NADH','FAD','UserChan'};

ValTypeLabs = {'abs. value','abs. change','% change'};
if BGShfBool ParLabs = {'shift','background','intensity','fraction engaged','\tau_1','\tau_2'}; end

close all

%% USER INPUT
OutFile = [''];
% OutFile = ['Analysis\OutName']; [a,b] = mkdir('Analysis'); % If sub-dir output
if size(mdpars,5)>1 mdpars(4,1,:,:,2,:,:,:) = 1-mdpars(4,1,:,:,2,:,:,:); end % FAD -> frac engaged
% Dims: [1=Param#, 2=mean/stderr(1,2), 3=time point, 4=position, 5=channel, 6=z-position, 7=embryo number, 8=segment]

% Exclusions
for i = 1:size(Excs,1)
    mdpars(:,:,:,Excs(i,1)+1,:,:,Excs(i,2),:)=nan;
end    

% Segments
Segs = {'joint','mito','cyto'};

% Choose which analysis to output
SubPlotsBool = 1; % All params in one fig w/ subplots
AvePlotsBool = 1; % If you want to plot the average curve in bold, enter '1'
IndivPlotsBool = 0; % Individual time plots with larger axes areas.

% Specify plot groupings and marker formatting
% col = [b;o];
col = colorblind; LineW = 1;
MarkS = 7;
LineTrans = 0.5;

% Specify time range and z-range. Leave empty to average all elements present
Trng = 1:size(mdpars,3);
Zrng = 1:size(mdpars,6);

%% Labels

Chlabs = {'NADH', 'FAD', 'User'};
% Detect which channels were present (not just filled with nan's)
ChTstmps = squeeze(nanmean(nanmean(nanmean(nanmean(mdpars(7,1,:,:,:,:,:),3),4),6),7));
chnums = find(~isnan(ChTstmps))';

% Derive number of exponents
npar = size(mdpars,1)-1; nexpo = floor((npar-2)/2); 
% Floor is because there is an extra fit param for synth IRFs.

% Make labels cells. IMPORTANT. Put them in the same order as how they are
% stored in the elements of MDpar. Then use pord to decide what order you
% want to plot them in.
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

%% AVERAGING: mainly to average across Z to get single value for each time point
% AVERAGES CALCULATED WITH WEIGHTING FACTORS = 1/SIGMA^2
% http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
% mdpars dims: [1=Param#, 2=mean/stderr(1,2), 3=time point, 4=position, 5=channel, 6=z-position, 7=embryo number, 8=segment]
% Fixed pars have 0 std err, causing NaNs in means. Just set those to 1.
mdpars(find(mdpars==0))=1;
if size(mdpars,6)==1
    ParZaved = mdpars;
else
    % Std error - from variance between z values.
    ParZaved = wtnanmean(mdpars(:,1,:,:,:,:,:,:),1./mdpars(:,2,:,:,:,:,:,:),6);
    ParZaved(:,2,:,:,:,:,:,:) = sqrt(wtnanvar(mdpars(:,1,:,:,:,:,:,:),1./mdpars(:,2,:,:,:,:,:,:),6))/sqrt(size(mdpars,6));
end

% % In case we want to look at Z dependence (we don't, generally)
if size(mdpars,3)==1
    ParTaved = mdpars;
else
    % Std error - from variance between z values.
    ParTaved = wtnanmean(mdpars(:,1,:,:,:,:,:,:),1./mdpars(:,2,:,:,:,:,:,:),3);
    ParTaved(:,2,:,:,:,:,:,:) = sqrt(wtnanvar(mdpars(:,1,:,:,:,:,:,:),1./mdpars(:,2,:,:,:,:,:,:),3))/sqrt(size(mdpars,3));
end

% Optional BG scaling operation. Generally not used, but sometimes useful check.
if BGsclBool
    load([acqpath 'BG_vals.mat'])
    % BGvals matrix dims: [ValType, framenum, posnum, ch];
    %    -> ValTypes: 1=intensity, 2=std, 3=NumOfPix, 4=timestamp
    
    load([acqpath 'multiD_indices.mat'])
    BadFrs = unique([nameinds{[nameinds{:,7}]==-1,6}]);
    ParZaved(:,:,BadFrs,:,:,:,:,:)=[];
    BGvals(:,BadFrs,:,:)=[];
    sz = size(ParZaved); 
    % Convert std to std err, reshape arrays
    BGstds = shiftdim(BGvals(2,:,:,:)./sqrt(BGvals(3,:,:,:)),1); % BGstds./sqrt(BGnum_pixels);
    BGstds = permute(repmat(BGstds,[1 1 1 1 sz(2) sz(6) sz(7) sz(8)]),[4 5 1 2 3 6 7 8]);
    BGirrs = shiftdim(BGvals(1,:,:,:),1);
    BGirrs = permute(repmat(BGirrs,[1 1 1 1 sz(2) sz(6) sz(7) sz(8)]),[4 5 1 2 3 6 7 8]);
    SclFct = nan(sz); SclFct(2:end,:,:,:,:,:,:,:)=[];
    SclFct(1,1,:,:,:,:,:,:) = BGirrs(:,1,:,:,:,:,:,:,:)./BGirrs(:,1,1,:,:,:,:,:,:);
    SclFct(1,2,:,:,:,:,:,:) = BGstds(:,1,:,:,:,:,:,:,:)./BGirrs(:,1,1,:,:,:,:,:,:); 
    % Scale intensities
    ParZaved(6,1,:,:,:,:,:,:) = ParZaved(6,1,:,:,:,:,:,:)./SclFct(1,1,:,:,:,1,:,:);
    ParZaved(6,2,:,:,:,:,:,:) = ...
    abs(ParZaved(6,1,:,:,:,:,:,:)).*sqrt( (ParZaved(6,2,:,:,:,:,:,:)./ParZaved(6,1,:,:,:,:,:,:)).^2 + (SclFct(1,2,:,:,:,:,:,:)./SclFct(1,1,:,:,:,:,:,:)).^2);
end

%% Time-sync'ing
% Perform optional time synchronization. Also convert to hours and start
% at 't=0' (not absolute time)
ParsTsnc = ParZaved;
% for rw = 1:size(ParsT,2)
%     for ch = 1:size(ParsTsnc,4)
%         for sg = 1:size(ParsTsnc,5)
%             if TmFrOffsetColm>-1
%                 TmFrOffset = ml{rw,TmFrOffsetColm};
%                 ParsTsnc(7,rw,:,ch,sg) = (ParsTsnc(7,rw,:,ch,sg)-ParsTsnc(7,rw,TmFrOffset,ch,sg))*86400/60/60;
%                 % Note, start time is 1st channel. FAD measurements come
%                 % slightly after NADH.
%             else % if no sync specified, just start each traj at 0
%                 ParsTsnc(7,rw,:,ch,sg) = (ParsTsnc(7,rw,:,ch,sg)-ParsTsnc(7,rw,1,ch,sg))*86400/60/60;
%             end
%         end
%     end
% end

%% Calculate average time curves

% Simplified from FLIM_mdpars_TimePlots, which interpolates. This routine
% averages masks within each position. Therefore no interpolation in
% necessary because all the masks within a position are imaged
% simultaneously.

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
% Dims of ParsTsnc: same as mdpars, but longer t-dim
% AllTs = ParsTsnc(7,1,:,:,:,:,:,:);
% T = min(AllTs(:)):(max(AllTs(:))-min(AllTs(:)))/1000:max(AllTs(:));
% if ~exist('TDispRan')|isempty(TDispRan) TDispRan = [min(T) max(T)]; end
% T(T<TDispRan(1)|T>TDispRan(2))=[]; % overwritten below
% sz = size(ParsTsnc);
% % Interp params
% ParsItrp = nan(sz(1),sz(2),size(T,2),sz(4),sz(5),sz(6),sz(7),sz(8));
% for p = 1:npar
%     for msk = 1:size(ParsTsnc,7)
%         for ch = 1:size(ParsTsnc,5)
%             for sg = 1:size(ParsTsnc,8)
%                 % Interp params
%                 x = squeeze(ParsTsnc(7,1,:,1,ch,1,msk,sg));
%                 y = squeeze(ParsTsnc(p,1,:,1,ch,1,msk,sg)); ystd = squeeze(ParsTsnc(p,2,:,1,ch,1,msk,sg));
%                 nanind = isnan(x) | isnan(y);
%                 x(nanind)=[]; y(nanind)=[]; ystd(nanind)=[];
%                 % Fill in dat colums with group values
%                 if ~isempty(x) % sometimes a segment is missing, e.g. no mito pixels for dim embryo
%                     ParsItrp(7,1,:,1,ch,1,msk,sg) = T;
%                     ParsItrp(p,1,:,1,ch,1,msk,sg) = interp1(x,y,T);
%                     % Also carry weights over to interpolation at discrete
%                     % time points for later plotting error bars.
%                     ParsItrp(p,2,:,1,ch,1,msk,sg) = interp1(x,ystd,T);
%                 end
%                 %                 % Test
%                 %                 plot(squeeze(ParsItrp(7,rw,:,ch,sg)),squeeze(ParsItrp(p,rw,:,ch,sg)));hold on;
%                 %                 plot(x,y,'.r');
%             end
%         end
%     end
% end

% Simplification, discrete positions rather than interp:
ParsItrp = ParsTsnc;


%% Load different calculates (plot types) into array to loop over
% Have multiple types of plots in this analysis:
% -> ValType: 1 = Average absolute values, 2 = Average absolute change,
%        3 = Average percentage change 
% CALCULATIONS
t0 = 1; % Frame of time point to divide time course by.
ParsCalc = nan([size(ParsItrp) 3]);
ParsCalc(:,:,:,:,:,:,:,:,1) = ParsItrp; % Abs

% Abs change
ParsCalc(:,:,:,:,:,:,:,:,2) = ParsItrp;
X = ParsItrp(1:6,1,:,:,:,:,:,:); dX = ParsItrp(1:6,2,:,:,:,:,:,:);
Y = ParsItrp(1:6,1,t0,:,:,:,:,:); dY = ParsItrp(1:6,2,t0,:,:,:,:,:);
F = X-Y; % Function - Relative change
ParsCalc(1:6,1,:,:,:,:,:,:,2) = F;
ParsCalc(1:6,2,:,:,:,:,:,:,2) = sqrt( dX.^2 + dY.^2 ); % Prop uncert

% Relative change
ParsCalc(:,:,:,:,:,:,:,:,3) = ParsItrp;
X = ParsItrp(1:6,1,:,:,:,:,:,:); dX = ParsItrp(1:6,2,:,:,:,:,:,:);
Y = ParsItrp(1:6,1,t0,:,:,:,:,:); dY = ParsItrp(1:6,2,t0,:,:,:,:,:);
F = X./Y; % Function - Relative change
ParsCalc(1:6,1,:,:,:,:,:,:,3) = F-1; 
ParsCalc(1:6,2,:,:,:,:,:,:,3) = abs(F).*sqrt( (dX./X).^2 + (dY./Y).^2); % Prop uncert


sz = size(ParsItrp);
GrpAves = nan([sz(1) sz(2) sz(3) sz(4) sz(5) sz(6) 1 sz(8) 3]);
GrpAves(1:6,1,:,:,:,:,:,:,:) = wtnanmean(ParsCalc(1:6,1,:,:,:,:,:,:,:),1./ParsCalc(1:6,2,:,:,:,:,:,:,:),7);
% GrpAves(1,:,:,:,:,:,:) = nanmean(ParsCalc(:,1,:,:,:,:,:,:),3); % Optional test against straight mean calc
GrpAves(1:6,2,:,:,:,:,:,:,:) = sqrt(wtnanvar(ParsCalc(1:6,1,:,:,:,:,:,:,:),1./ParsCalc(1:6,2,:,:,:,:,:,:,:),7))./sqrt(sz(7)); % Std errs
GrpAves(7,:,:,:,:,:,:,:,:,:) = ParsCalc(7,:,:,:,:,:,1,:,:); % Times shouldn't get averaged, just transferred

% % Test
% plot(squeeze(GrpAves(7,1,:,ch,sg)),squeeze(GrpAves(4,1,:,ch,sg)));hold on;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% TIME PLOTS

% Determine if there are multiple time points:
ts = unique(nanmean(nanmean(nanmean(ParsCalc(7,1,:,:,:,:,1,1),4),5),6)); ts(isnan(ts))=[];

if length(ts)>1
% Plot all 8 params on single plot with Subplots
for vt = ValInds % value type -> different calculations derived in ParsCalc
    if SubPlotsBool
        % Set figure sizes
        MidW = .15; CapS = 20; Crop = 0.7; pos = [51 113 1300 684];
        Spread = 13; FontS = 16;
        if ~BGShfBool pord = [6 4 3 5]; pos = [51 113 1600 684]; end
        close all;
        
        % Define a time gap to plot spaced out, interpolated error bars
        %     Tgap = round(1001/size(ParsTsnc,3))*1.7; % 1.7 factor to make them look nicely spaced
        
        for posnum = 1:size(mdpars,4)
            for seg = 1:size(mdpars,8)
                subind = 1;
                for p = 1:length(pord)
                    for ch = Chind
                        subplot(length(chnums),length(pord),p+length(pord)*(ch-1))
                        %             dat = nan(size(mdpars,7),length(T));
                        % Fill in dat colums with group values
                        dat = squeeze(ParsCalc(pord(p),1,:,posnum,ch,1,:,seg,vt));
                        daterr = squeeze(ParsCalc(pord(p),2,:,posnum,ch,1,:,seg,vt));
                        T = squeeze(ParsCalc(end,1,:,posnum,ch,1,:,seg,vt));
                        T = (T - min(T)).*86400./timeunit;
                        TDispRan = [min(T(:)) max(T(:))];
                        % Dims of ParsCalc: [par, ListRow, time(synth), channel, seg]
                        if vt==3 daterr = daterr*100; dat = dat*100; end % Convert to percent
                        if ~isempty(dat(~isnan(dat)))
                        plts = plot(T,dat,'color',[col(1,:) LineTrans],'linewidth',LineW); hold on;
                        pltsE = errorbar(gca,T,dat,daterr,'color',[col(1,:) LineTrans],'linestyle','none'); hold on;
%                         for pl = 1:length(plts)
%                             plts(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
%                             plts(pl).Color = [col(pl,:) LineTrans];
%                             pltsE(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
%                             pltsE(pl).Color = [col(pl,:) LineTrans];
%                             PltLbs{pl} = ['m' num2str(pl)];
%                         end
                        
                        for pl = 1:size(dat,2)
                            x = T(:,pl); y = dat(:,pl); yerr = daterr(:,pl); 
                            rm = isnan(x); x(rm)=[]; y(rm)=[]; yerr(rm)=[];
                            if ~isempty(x)
                            plts(pl) = plot(x,y,'color',[col(1,:) LineTrans],'linewidth',LineW); hold on;
                            pltsE(pl) = errorbar(x,y,yerr,'color',[col(1,:) LineTrans],'linestyle','none'); hold on;
                            plts(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                            plts(pl).Color = [col(pl,:) LineTrans];
                            pltsE(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                            pltsE(pl).Color = [col(pl,:) LineTrans];
                            PltLbs{pl} = ['m' num2str(pl)];
                            end
                        end
                        
                        h(1) = plts(1);
                        GrpAveX = nanmean(T,2); GrpAveY = squeeze(GrpAves(pord(p),1,:,posnum,ch,1,1,seg,vt));
                        GrpAveErr = squeeze(GrpAves(pord(p),2,:,posnum,ch,1,1,seg,vt));
                        rm = isnan(GrpAveX); GrpAveX(rm)=[]; GrpAveY(rm)=[]; GrpAveErr(rm)=[];
                        %             if AvePlotsBool
                        if vt==3 GrpAveErr = GrpAveErr*100; GrpAveY = GrpAveY*100; end % Convert to percent
                        errorbar(gca,GrpAveX,GrpAveY,GrpAveErr,'color','k','LineWidth',1.5);
                        %                 h(1) = plot(GrpAveX,GrpAveY,'color',col(1,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_Grp' num2str(Grps{1})]);
                        %             end
                        
                        TightWSpace(gca,.05)
                        xlim(TDispRan)
                        set(gca,'fontsize',FontS); set(gcf,'color','w')
                        end
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
                %legend(plts,PltLbs,'location','best','fontsize',8);
                xlabel(tlab,'position',[.32 -.142],'fontsize',FontS)
                legend(plts,PltLbs,'position',[0.92521    0.0042981     0.073125       0.2288]);
                flab = ValTypeLabs{vt}; flab(strfind(flab,'.'))=[]; flab(strfind(flab,' '))=[];
                if BGsclBool flab = [flab '_BGscl']; end
                if BGShfBool flab = [flab '_BGShft']; end
                saveas(gcf,[acqpath 'Tsubplots_Pos' num2str(posnum-1) '_' flab '_' Segs{seg} '_' MDlab '.png']);
                saveas(gcf,[acqpath 'Tsubplots_Pos' num2str(posnum-1) '_' flab '_' Segs{seg} '_' MDlab '.fig']);
            end
        end
    end
end
end

%% Z-PLOTS

% Determine if there are multiple z-planes. If so, plot as a function of
% frame. Also create an additional x-axis that shows the physical z
% position (taken from 'AcqSettings.txt')

if size(mdpars,6)>3 % Generally time courses only have 3 planes, max
AcqSet = ReadAcqSettingsTxt([acqpath 'AcqSettings.txt']); Zcell = AcqSet(end,:);
for i=1:length(Zcell)
    if ischar(Zcell{i})
        if ~isempty(Zcell{i}) Zs(i) = str2num(Zcell{i}); end
    else
        if ~isempty(Zcell{i}) Zs(i) = Zcell{i}; end
    end
end
%     % test two time points, two masks
%     md2 = mdpars; md2(:,1,:,:,:,:,:,:)=md2(:,1,:,:,:,:,:,:)+1; mdpars = cat(3,mdpars,md2); 
%     md2 = mdpars; md2(:,1,:,:,:,:,:,:)=md2(:,1,:,:,:,:,:,:)-1; mdpars = cat(7,mdpars,md2); 
Xlbs = {'Frame #','Z (\mum)'};Xflbs = {'XaxFr','XaxZum'};

for xlb = 1%:2 % 1 plots with 'frames' x-axis. Set to 2 to plot vs microns
for vt = ValInds % value type -> different calculations derived in ParsCalc
    if SubPlotsBool
        % Set figure sizes
        MidW = .15; CapS = 20; Crop = 0.7; pos = [51 113 1300 684];
        Spread = 13; FontS = 16;
        if ~BGShfBool pord = [6 4 3 5]; pos = [51 113 1600 684]; end
        close all;
        
        % Define a time gap to plot spaced out, interpolated error bars
        %     Tgap = round(1001/size(ParsTsnc,3))*1.7; % 1.7 factor to make them look nicely spaced
        
        for posnum = 1:size(mdpars,4)
            for seg = 1:size(mdpars,8)
                subind = 1;
                for p = 1:length(pord)
                    for ch = Chind
                        subplot(length(chnums),length(pord),p+length(pord)*(ch-1))
                        %             dat = nan(size(mdpars,7),length(T));
                        % Fill in dat colums with group values
                        dat = permute(mdpars(pord(p),1,:,posnum,ch,:,:,seg,vt),[6 3 7 1 2 4 5 8]);
                        dat = reshape(dat,size(mdpars,3)*size(mdpars,6),size(mdpars,7)); 
                        % Frames go down rows, columns are separate masks, if present
                        daterr = permute(mdpars(pord(p),2,:,posnum,ch,:,:,seg,vt),[6 3 7 1 2 4 5 8]);
                        daterr = reshape(daterr,size(mdpars,3)*size(mdpars,6),size(mdpars,7));
%                         % Time stamps, but don't use these. Instead use frames and z-positions
%                         T = permute(mdpars(7,1,:,posnum,ch,:,:,seg,vt),[6 3 7 1 2 4 5 8]);
%                         T = reshape(T,size(mdpars,3)*size(mdpars,6),size(mdpars,7)); 
%                         T = (T - min(T)).*86400./timeunit;
%                         TDispRan = [min(T(:)) max(T(:))];
                        if xlb==1
                            Xarr = repmat((1:size(dat,1))',1,size(mdpars,7)); % Frames
                        else
                            Xarr = repmat(Zs',size(mdpars,3),size(mdpars,7)); % z-plane in microns
                        end
                        % Dims of ParsCalc: [par, ListRow, time(synth), channel, seg]
                        if vt==3 daterr = daterr*100; dat = dat*100; end % Convert to percent
                        if ~isempty(dat(~isnan(dat)))
                        plts = plot(Xarr,dat,'color',[col(1,:) LineTrans],'linewidth',LineW); hold on;
                        pltsE = errorbar(gca,Xarr,dat,daterr,'color',[col(1,:) LineTrans],'linestyle','none'); hold on;
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
                            plts(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                            plts(pl).Color = [col(pl,:) LineTrans];
                            pltsE(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                            pltsE(pl).Color = [col(pl,:) LineTrans];
                            PltLbs{pl} = ['m' num2str(pl)];
                            end
                        end
                        
                        h(1) = plts(1);
                        % Don't need average plots for these.
%                         GrpAveX = nanmean(Xarr,2); GrpAveY = squeeze(GrpAves(pord(p),1,:,posnum,ch,1,1,seg,vt));
%                         GrpAveErr = squeeze(GrpAves(pord(p),2,:,posnum,ch,1,1,seg,vt));
%                         rm = isnan(GrpAveX); GrpAveX(rm)=[]; GrpAveY(rm)=[]; GrpAveErr(rm)=[];
%                         %             if AvePlotsBool
%                         if vt==3 GrpAveErr = GrpAveErr*100; GrpAveY = GrpAveY*100; end % Convert to percent
%                         errorbar(gca,GrpAveX,GrpAveY,GrpAveErr,'color','k','LineWidth',1.5);
                        %                 h(1) = plot(GrpAveX,GrpAveY,'color',col(1,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_Grp' num2str(Grps{1})]);
                        %             end
                        
                        TightWSpace(gca,.05)
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
                legend(plts,PltLbs,'location','best','fontsize',8);
                xlabel(Xlbs{xlb},'position',[.32 -.142],'fontsize',FontS)
                legend(plts,PltLbs,'position',[0.92521    0.0042981     0.073125       0.2288]);
                flab = ValTypeLabs{vt}; flab(strfind(flab,'.'))=[]; flab(strfind(flab,' '))=[];
                if BGsclBool flab = [flab '_BGscl']; end
                if BGShfBool flab = [flab '_BGShft']; end
                saveas(gcf,[acqpath 'Zsubplots_Pos' num2str(posnum-1) '_' flab '_' Segs{seg} '_' MDlab '_' Xflbs{xlb} '.png']);
                saveas(gcf,[acqpath 'Zsubplots_Pos' num2str(posnum-1) '_' flab '_' Segs{seg} '_' MDlab '_' Xflbs{xlb} '.fig']);
            end
        end
    end
end
end
end
end
