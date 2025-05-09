% Running 'FLIM_master_list_and_mdpars' creates a master list and multiD
% parameter matrix. Copy this script and those list/MDpar files to the
% anlaysis folder you want, and input the file names below. Then this
% script will make plots and do analysis on the data.

% This function is designed explicitly for concatenating datasets along the
% time dimension. Sometimes acquisitions get interrupted and have to be
% restarted, splitting the time courses into separate folders. In the
% Master List, all masks with the same 'Data Set' number will be
% concatenated together.

% NOTE: All of these operations treat each sample (embryo/oocytes) as a
% single data point. That means averaging over time points and z-planes,
% but you can specify what Trng and Zrng you want below.

clear all
cd(pwd) % Place a copy of this script in the folder where you want plots to export.
close all
load('Z:\Lab\Tim\MATLAB TIM\Metabolic Imaging\Plotting and analysis\colorblind_colormap.mat')

%% USER INPUT
OutFile = ['CatTrajs'];
% OutFile = ['Analysis\OutName']; [a,b] = mkdir('Analysis'); % If sub-dir output
% Load and concatenate lists
listpaths{1} = 'MasterList_1stIVM.xls';
% listpaths{2} = 'MasterList_1.xls';
parpaths{1} = 'MasterList_1stIVM_pars.mat'; 
% parpaths{2} = 'MasterList_2.mat';
[ml,MDpars,txt] = GroupMasterListPars(listpaths,parpaths);
% Convert FAD long fraction to frac 'engaged', inverse
MDpars(4,1,:,:,2,:,:) = 1-MDpars(4,1,:,:,2,:,:);

% Choose segment
seg = 1; Segs = {'joint','mito','cyto'};

% Choose which analysis to output
SubPlotsBool = 1; % All params in one fig w/ subplots
IndivPlotsBool = 0; % Individual time plots with larger axes areas.

% Specify plot groupings and marker formatting
Grps = {[1 2 3 5],[4 6]}; % Group sets in cells to plot together, e.g. Grps = {[2 4],[1 3]}
markers = {'o','o','o','o','o','o','o','o'};
r=[1 0 0]; g=[0 1 0]; b=[0 0 1]; bl=[117 255 255]./255; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1];
% col = [b;o];
col = colorblind;
LineW = 1;
MarkS = 7;
LineTrans = 0.4;

% Set time range. 
% 1) [] Default  - plot full time range, including far time points that only
%               have few or only one trajectory. 
%               NOTE: average curves at these extremes show large
%               fluctuations due to small numbers of trajectories.
%               Typically, you want to crop at some width where you have
%               plenty of trajectories and the average is well-behaved.
% 2) [Tmin,Tmax] (array of doubles) - set custom time range
TDispRan = [];
% Note, for dev curves, [-5 65] good for 1st div, [-22 40] good for 2nd
% div, and [-65 5] for blastocoel

% Time frame offset. Enter -1 if you don't want to shift time courses.
% 6 is the 1st division, 7 is 2nd division, 8 is blastocoel
TmFrOffsetColm = -1;

% If you want to plot the average curve in bold, enter '1'
AvePlotsBool = 1;

% Specify time range and z-range. Leave empty to average all elements present
Trng = []; if isempty(Trng) Trng = 1:size(MDpars,4); end
Zrng = []; if isempty(Zrng) Zrng = 1:size(MDpars,6); end

%% Auto-setup groups and labels
% Get the group numbers in an array for easy reference. 
% Normally, get Grp inds (row indices) for each group. But now Data Set
% numbers specify which trajectories to concatenate. Thus, each set should
% only have 1 element (trajectory) after concatenating, and each 'group'
% should have only the indices of the individual trajectories specified by
% user. 
% To concatenate, define a different cell of indices, 'Usets' for 'unique
% data sets detected in spreadsheet.
SetsArr = [ml{:,2}]; Usets = unique(SetsArr); 
Usets(isnan(Usets))=[]; Usets = sort(Usets);
for dset = 1:length(Usets)
    Usetinds{dset} = find(SetsArr==Usets(dset));
    SLs(dset) = length(Usetinds{dset}); % Array lengths of each group
end

% 'row' indices of concatenated MDpars array will correspond to data sets.
% Therefor 'Grpinds' will be equal to the the 'Grps' cell specified by the
% user.
Grpinds = Grps;
% Get group labels from first of grouped sets
for grp = 1:length(Grps)
    ind = find(SetsArr==Grpinds{grp}(1));
    GrpLabs{grp} = ml{ind(1),3};
    Ls(grp) = length(Grpinds{grp}); % Array lengths of each group
end
Chlabs = {'NADH', 'FAD', 'User'};

npar = size(MDpars,1)-1; nexpo = (npar-2)/2;
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
% MDpars dims: [Param#, mean/std(1,2), ListRow, time point, channel, z-position, segment]
% Fixed pars have 0 std err, causing NaNs in means. Just set those to 1.
MDpars(find(MDpars==0))=1;
ParZaved = wtnanmean(MDpars(:,1,:,Trng,:,Zrng,:),1./MDpars(:,2,:,Trng,:,Zrng,:),6);
ParZaved(:,2,:,:,:,:,:) = sqrt(1./nansum(1./MDpars(:,2,:,Trng,:,Zrng,:).^2,6));
ParsT = squeeze(ParZaved(:,1,:,:,:,:,:)); % Take par vals only
% Dims of ParsT: [par, ListRow, time, channel, seg].

% In case we want to look at Z dependence
ParTaved = wtnanmean(MDpars(:,1,:,Trng,:,Zrng,:),1./MDpars(:,2,:,Trng,:,Zrng,:),4);
ParTaved(:,2,:,:,:,:,:) = sqrt(1./nansum(1./MDpars(:,2,:,Trng,:,Zrng,:).^2,4)); % Propagated error
ParsZ = squeeze(ParTaved);

%% Concatenate set param arrays along time dim (for each specified group)
for dset = 1:length(Usetinds)
    Pars2Cat = ParZaved(:,:,Usetinds{dset},:,:,:,:);
    CatPars = Pars2Cat(:,:,1,~isnan(Pars2Cat(3,1,1,:,1,1,1)),:,:,:);
    sz(1) = size(CatPars,4);
    for i = 2:size(Pars2Cat,3)
        nanind = squeeze(~isnan(Pars2Cat(3,1,i,:,1,1,1)));
        sz(i)=size(Pars2Cat(:,:,i,nanind,:,:,:),4);
        CatPars = cat(4,CatPars,Pars2Cat(:,:,i,nanind,:,:,:));
    end
    CatLs(dset) = size(CatPars,4);
    CatCell{dset} = CatPars;
end
sz = size(ParZaved);
MDcat = nan([sz(1) sz(2) length(Usetinds) max(CatLs) sz(5) sz(6) sz(7)]);
% Make each concatenated group into a single element of a new MD matrix
for dset = 1:length(Usetinds)
    MDcat(:,:,dset,1:CatLs(dset),:,:,:) = CatCell{dset};
end

% % Test
% size(MDcat)
% plot(squeeze(MDcat(7,1,11,:,1,1,1)),squeeze(MDcat(6,1,11,:,1,1,1)))



%% Time-sync'ing - NOT IN THIS VERSION YET
% Perform optional time synchronization. Also convert to hours and start
% at 't=0' (not absolute time)
ParsTsnc = squeeze(MDcat(:,1,:,:,:,:,:));
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
AllTs = ParsTsnc(7,:,:,:,:);
T = min(AllTs(:)):(max(AllTs(:))-min(AllTs(:)))/1000:max(AllTs(:));
% T hasn't been converted to hours yet because that was previously done in
% T-sync. Do it here instead.
T = (T-T(1))*86400/60/60;
sz = size(ParsTsnc);
% Interp params
ParsItrp = nan(sz(1),sz(2),size(T,2),sz(4),sz(5));
for p = 1:npar
    for rw = 1:size(ParsTsnc,2)
        for ch = 1:size(ParsTsnc,4)
            for sg = 1:size(ParsTsnc,5)
                % Interp params
                x = squeeze(ParsTsnc(7,rw,:,ch,sg));
                y = squeeze(ParsTsnc(p,rw,:,ch,sg));
                nanind = isnan(x) | isnan(y);
                x(nanind)=[]; y(nanind)=[];
                x = (x-x(1))*86400/60/60;
                % Fill in dat colums with group values
                ParsItrp(7,rw,:,ch,sg) = T;
                ParsItrp(p,rw,:,ch,sg) = interp1(x,y,T);
                %                 % Test
                %                 plot(squeeze(ParsItrp(7,rw,:,ch,sg)),squeeze(ParsItrp(p,rw,:,ch,sg)));hold on;
                %                 plot(x,y,'.r');
            end
        end
    end
end
% Get average curves for groups
GrpAves = nan(sz(1),length(Grpinds),size(T,2),sz(4),sz(5));
for grp = 1:length(Grpinds)
    GrpAves(:,grp,:,:,:) = nanmean(ParsItrp(:,Grpinds{grp},:,:,:),2);
end
% % Test
% plot(squeeze(GrpAves(7,1,:,ch,sg)),squeeze(GrpAves(4,1,:,ch,sg)));hold on;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Plot all 8 params on single plot with Subplots
if SubPlotsBool
    % Set figure sizes
    MidW = .15; CapS = 20; Crop = 0.7; pos = [800 100 408 485];
    Spread = 13; FontS = 16;
    pord = [6 4 3 5];
%     pord = [1 2 6 4 3 5]; % If you want to include shift and BG
    if ~exist('TDispRan')|isempty(TDispRan) TDispRan = [min(T) max(T)]; end
    close all;
    for p = 1:length(pord)
        for ch = 1:2
            subplot(2,length(pord),p+length(pord)*(ch-1))
            for grp = 1:length(Grpinds)
                dat = nan(Ls(grp),length(T));
                % Fill in dat colums with group values
                dat(1:length(Grpinds{grp}),:) = ParsItrp(pord(p),Grpinds{grp},:,ch,seg);
                % Dims of ParsItrp: [par, ListRow, time(synth), channel, seg]
                Lab = Labs{p}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
                
                plot(T,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW,'DisplayName',[Chlabs{ch} '_' GrpLabs{grp}]); hold on;
                GrpAveX = squeeze(GrpAves(7,grp,:,ch,seg)); GrpAveY = squeeze(GrpAves(pord(p),grp,:,ch,seg));
                h(grp) = plot(GrpAveX,GrpAveY,'color',col(grp,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_' Lab]);
            end
            TightWSpace(gca,.05)
            xlim(TDispRan)
            set(gca,'fontsize',FontS); set(gcf,'color','w')
            xlabel('Time (h)');
            ylabel([Chlabs{ch} ' ' Labs{pord(p)}],'fontsize',FontS*1.1)
        end
    end
    set(gcf,'position',[51 113 1300 684])
    ResizeSubplots(gcf,[.05 .05 .02 .04],[.05 .07])
    legend(h,GrpLabs,'location','best');
    saveas(gcf,[OutFile '_Tsubplots_' Segs{seg} '.png']);
end


%% Individual Time plots
if IndivPlotsBool 
    % Set figure sizes
    MidW = .15; CapS = 20; Crop = 0.7; pos = [420   278   560   480];
    Spread = 13; FontS = 16;
    
    if ~exist('TDispRan') TDispRan = [min(T) max(T)]; end
    close all
    
    for p = 1:npar
        for ch = 1:2
            close all;
            for grp = 1:length(Grpinds)
                dat = nan(Ls(grp),length(T));
                % Fill in dat colums with group values
                dat(1:length(Grpinds{grp}),:) = ParsItrp(p,Grpinds{grp},:,ch,seg);
                % Dims of ParsItrp: [par, ListRow, time(synth), channel, seg]
                Lab = Labs{p}; if strcmp(Lab(1),'\') Lab = Lab(2:end); end
                
                
                plot(T,dat,'color',[col(grp,:) LineTrans],'linewidth',LineW,'DisplayName',[Chlabs{ch} '_' GrpLabs{grp}]); hold on;
                GrpAveX = squeeze(GrpAves(7,grp,:,ch,seg)); GrpAveY = squeeze(GrpAves(p,grp,:,ch,seg));
                plot(GrpAveX,GrpAveY,'color',col(grp,:),'linewidth',3,'DisplayName',[Chlabs{ch} '_' Lab]);
            end
            TightWSpace(gca,.05)
            xlim(TDispRan)
            set(gca,'fontsize',FontS); set(gcf,'color','w')
            xlabel('Time (h)');
            ylabel([Chlabs{ch} ' ' Labs{p}],'fontsize',FontS*1.1)
            
            set(gcf,'position',pos)
            axpos = get(gca,'position'); set(gca,'position',axpos+[0 0 0.06 0])
            Lab(Lab==' ')=[]; % Remove spaces
            saveas(gcf,[OutFile '_Tplot_' Chlabs{ch} '_' Lab '_' Segs{seg} '.png']);
        end
    end
end