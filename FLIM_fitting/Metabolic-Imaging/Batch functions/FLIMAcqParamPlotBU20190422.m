function FLIMAcqParamPlot(acqpath,Excs)
% Plots params from a saved data set from mdpars (multi-d parameter matrix)
% Adapted from FLIM_mdpars_TimePlots to act on individual positions of
% acquisition. Quicker way to visualize acquisitions that doing plots for
% individual masks.
% INPUTS:
% -acqpath: duh
% -Excs: exclusions. nx2 array of pos-mask pairs. Can rerun this after
%   initial analysis to exclude bad masks.

% clear all;
% acqpath = 'C:\Dropbox\data\s2_a1\';
% Excs = [0 1];

if acqpath(end)~='\' acqpath = [acqpath '\']; end
if ~exist('Excs') Excs = []; end
load('C:\Google Drive\MATLAB PROGRAMS\Metabolic Imaging\Plotting and analysis\colorblind_colormap.mat')
colorblind = [colorblind; rand(10,3)]; % Add rands in case of more than 12 masks in one position (rare)
colorblind(5,:)=[]; % Delete black, that's used for averages
% Load mdpars
load([acqpath 'multiD_pars.mat']);

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

close all

%% USER INPUT
OutFile = [''];
% OutFile = ['Analysis\OutName']; [a,b] = mkdir('Analysis'); % If sub-dir output
mdpars(4,1,:,:,2,:,:,:) = 1-mdpars(4,1,:,:,2,:,:,:); % FAD -> frac engaged
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
BGShfBool = 0; % Whether to display BG and shift on subplots

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

npar = size(mdpars,1)-1; nexpo = (npar-2)/2;
% Make labels cells. IMPORTANT. Put them in the same order as how they are
% stored in the elements of MDpar. Then use pord to decide what order you
% want to plot them in.
if nexpo==2
    % Short labels
    Labs = {'Shift','A(1-bg)','tau1','frac2','tau2','int'};
    % Labels for Seli paper
    Labs = {'Shift','A(1-bg)','\tau_1 (ns)','Fraction engaged','\tau_2 (ns)','Intensity (au)'};
    pord = [ 1 2 6 3 4 5];
elseif nexpo==3
    Labs = {'Shift','A(1-bg)','tau1','frac2','tau2','frac3','tau3','int'};
    pord = [1 2 8 3 4 5 6 7];
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
    ParZaved(:,2,:,:,:,:,:) = sqrt(wtnanvar(MDpars(:,1,:,Trng,:,Zrng,:),1./MDpars(:,2,:,Trng,:,Zrng,:),6))/sqrt(size(MDpars,6));
    
end
% Previously just:
%     ParZaved = wtnanmean(mdpars(:,1,:,Trng,:,Zrng,:),1./mdpars(:,2,:,Trng,:,Zrng,:),6);
%     ParZaved(:,2,:,:,:,:,:) = sqrt(1./nansum(1./mdpars(:,2,:,Trng,:,Zrng,:).^2,6));

% Dims of ParsT: [par, ListRow, time, channel, seg].

% In case we want to look at Z dependence (we don't, generally)
ParTaved = wtnanmean(mdpars(:,1,:,:,:,:,:),1./mdpars(:,2,:,:,:,:,:),3);
ParTaved(:,2,:,:,:,:,:) = sqrt(1./nansum(1./mdpars(:,2,:,:,:,:,:).^2,3)); % Propagated error
ParsZ = squeeze(ParTaved);

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
%        3 = Average relative change
t0 = find(abs(T)==min(abs(T)));
ParsCalc = nan([size(ParsItrp) 3]);
ParsCalc(:,:,:,:,:,:,:,1) = ParsItrp; % Abs
ParsCalc(:,:,:,:,:,:,:,2) = ParsItrp;
ParsCalc(1:6,1,:,:,:,:,:,2) = ParsItrp(1:6,1,:,:,:,:,:)-ParsItrp(1:6,1,:,t0,:,:,:); % Abs change
ParsCalc(1:6,2,:,:,:,:,:,2) = sqrt(ParsItrp(1:6,2,:,:,:,:,:).^2+ParsItrp(1:6,2,:,t0,:,:,:).^2); % prop of std error
ParsCalc(:,:,:,:,:,:,:,3) = ParsItrp;
ParsCalc(1:6,1,:,:,:,:,:,3) = ParsItrp(1:6,1,:,:,:,:,:)./ParsItrp(1:6,1,:,t0,:,:,:)-1; % Relative change
% I tried to do std err prop on the relative one using this way:
ParsCalc(1:6,2,:,:,:,:,:,3) = abs(ParsItrp(1:6,1,:,:,:,:,:)).*sqrt( (ParsItrp(1:6,2,:,:,:,:,:)./ParsItrp(1:6,1,:,:,:,:,:)).^2 + (ParsItrp(1:6,2,:,t0,:,:,:)./ParsItrp(1:6,1,:,t0,:,:,:)).^2); %
% But it kept doing something goofy and causing an error with the means
% Instead just divide the std err at tR by the initial mean as an approx
% ParsCalc(1:6,2,:,:,:,:,:,3) = ParsItrp(1:6,2,:,:,:,:,:)./abs(ParsItrp(1:6,1,:,t0,:,:,:));




% Get average curves for this position
sz = size(ParsItrp);
GrpAves = nan([sz(1) sz(2) sz(3) sz(4) sz(5) sz(6) 1 sz(8)]);
GrpAves(1:6,1,:,:,:,:,:,:) = wtnanmean(ParsItrp(1:6,1,:,:,:,:,:,:),1./ParsItrp(1:6,2,:,:,:,:,:,:),7);
% GrpAves(1,:,:,:,:,:,:) = nanmean(ParsItrp(:,1,:,:,:,:,:,:),3); % Optional test against straight mean calc
GrpAves(1:6,2,:,:,:,:,:,:) = sqrt(wtnanvar(ParsItrp(1:6,1,:,:,:,:,:,:),1./ParsItrp(1:6,2,:,:,:,:,:,:),7))/sqrt(size(ParsItrp,7)); % Std errs
GrpAves(7,1,:,:,:,:,:,:) = ParsItrp(7,1,:,:,:,:,1,:); % Times shouldn't get averaged, just transferred

% % Test
% plot(squeeze(GrpAves(7,1,:,ch,sg)),squeeze(GrpAves(4,1,:,ch,sg)));hold on;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Plot all 8 params on single plot with Subplots
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
            
            for p = 1:length(pord)
                for ch = 1:2
                    
                    subplot(2,length(pord),p+length(pord)*(ch-1))
                    %             dat = nan(size(mdpars,7),length(T));
                    % Fill in dat colums with group values
                    dat = squeeze(ParsItrp(pord(p),1,:,posnum,ch,1,:,seg));
                    daterr = squeeze(ParsItrp(pord(p),2,:,posnum,ch,1,:,seg));
                    T = squeeze(ParsItrp(7,1,:,posnum,ch,1,:,seg));
                    T = (T - min(T)).*86400./timeunit;
                    TDispRan = [min(T(:)) max(T(:))];
                    % Dims of ParsItrp: [par, ListRow, time(synth), channel, seg]
                    Lab = Labs{p}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
                    
                    plts = plot(T,dat,'color',[col(1,:) LineTrans],'linewidth',LineW); hold on;
                    pltsE = errorbar(gca,T,dat,daterr,'color',[col(1,:) LineTrans],'linestyle','none'); hold on;
                    for pl = 1:length(plts) 
                        plts(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)]; 
                        plts(pl).Color = [col(pl,:) LineTrans];
                        pltsE(pl).DisplayName = [Chlabs{ch} '_Pos' num2str(posnum-1) '_' num2str(pl)];
                        pltsE(pl).Color = [col(pl,:) LineTrans];
                        PltLbs{pl} = ['m' num2str(pl)];
                    end
                    h(1) = plts(1);
                    GrpAveX = nanmean(T,2); GrpAveY = squeeze(GrpAves(pord(p),1,:,posnum,ch,1,1,seg));
                    GrpAveErr = squeeze(GrpAves(pord(p),2,:,posnum,ch,1,1,seg));
                    %             if AvePlotsBool
                    errorbar(gca,GrpAveX,GrpAveY,GrpAveErr,'color','k','LineWidth',1.5);
                    %                 h(1) = plot(GrpAveX,GrpAveY,'color',col(1,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_Grp' num2str(Grps{1})]);
                    %             end
                    
                    TightWSpace(gca,.05)
                    xlim(TDispRan)
                    set(gca,'fontsize',FontS); set(gcf,'color','w')
                    xlabel(tlab);
                    ylabel([Chlabs{ch} ' ' Labs{pord(p)}],'fontsize',FontS*1.1)
                end
            end
            set(gcf,'position',pos)
            ResizeSubplots(gcf,[.04 .05 .02 .04],[.05 .07])
            legend(plts,PltLbs,'location','best','fontsize',8);
            saveas(gcf,[acqpath 'Tsubplots_Pos' num2str(posnum-1) '_' Segs{seg} '.png']);
            saveas(gcf,[acqpath 'Tsubplots_Pos' num2str(posnum-1) '_' Segs{seg} '.fig']);
        end
    end
end
    
%     %% Individual Time plots
%     if IndivPlotsBool
%         % Set figure sizes
%         MidW = .15; CapS = 20; Crop = 0.7; pos = [420   278   560   480];
%         Spread = 13; FontS = 16;
%         
%         if ~exist('TDispRan') TDispRan = [min(T) max(T)]; end
%         close all
%         
%         for p = 1:npar
%             for ch = 1:2
%                 close all;
%                 for grp = 1:length(Grpinds)
%                     dat = nan(Ls(grp),length(T));
%                     % Fill in dat colums with group values
%                     dat(1:length(Grpinds{grp}),:) = ParsItrp(p,Grpinds{grp},:,ch,seg);
%                     % Dims of ParsItrp: [par, ListRow, time(synth), channel, seg]
%                     Lab = Labs{p}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
%                     
%                     
%                     plts = plot(T,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW,'DisplayName',[Chlabs{ch} '_' GrpLabs{grp}]); hold on;
%                     h(grp) = plts(1);
%                     GrpAveX = squeeze(GrpAves(7,grp,:,ch,seg)); GrpAveY = squeeze(GrpAves(p,grp,:,ch,seg));
%                     if AvePlotsBool
%                         h(grp) = plot(GrpAveX,GrpAveY,'color',col(grp,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_' Lab]);
%                     end
%                 end
%                 TightWSpace(gca,.05)
%                 xlim(TDispRan)
%                 set(gca,'fontsize',FontS); set(gcf,'color','w')
%                 xlabel('Time (min)');
%                 ylabel([Chlabs{ch} ' ' Labs{p}],'fontsize',FontS*1.1)
%                 
%                 set(gcf,'position',pos)
%                 axpos = get(gca,'position'); set(gca,'position',axpos+[0 0 0.06 0])
%                 Lab(Lab==' ')=[]; % Remove spaces
%                 legend(h,GrpLabs,'location','best');
%                 saveas(gcf,[OutFile '_Tplot_' Chlabs{ch} '_' Lab '_' Segs{seg} '.png']);
%             end
%         end
%     end
    
    
