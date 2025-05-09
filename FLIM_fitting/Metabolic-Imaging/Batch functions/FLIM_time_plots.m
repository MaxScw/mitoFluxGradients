
% Running 'FLIM_batch_averages' calculates all the average params
% in a folder and outputs values to an excell sheet. Enter
% the row ranges and species names below and this script will create
% various plots of the parameters:
% - NAD Irr vs FAD Irr
% - NAD Fbound vs FAD Fbound
% - NAD Tau1 vs FAD Tau1
% - NAD Tau2 vs FAD Tau2
% - NAD Irr vs FAD Irr
clear all

path = pwd;
path = UpOneDir(path);
% load([path 'TimeArrs'])
suffix = 'Pos3MaskComp';
OutFile = [path '\Analysis\FLIM_time_plots' suffix];
AvePlotsBool = 1;
% Optional legend labels
% GroupNames = {'Blast5%','Dead cells5%','Whole Emb5%','Blast2%','Dead cells2%','Whole Emb2%'};
GroupNames = {'Blasts','Weak Blasts','Frag'};
%% USER INPUT
OutFile = ['Z:\Lab\Tim\2017-03-10 Racowsky - Human Embs 5 and 2 per O2\Analysis and Plots\Bothsets_samecolor_BlastsWeakFrag'];

% Data set 1
path = 'Z:\Lab\Tim\2017-03-10 Racowsky - Human Embs 5 and 2 per O2\';
file = [path 'ParamsAllSams_alltime.xls'];
Tfile = [path 'TimeArrs'];
[num,txt,res]= xlsread(file);
% Remove 1st row of headers
txt(1,:)=[]; res(1,:)=[];
dat = load(Tfile);
TimeArrs = dat.TimeArrs; TimeArrs_std = dat.TimeArrs_std;


% Data set 2 (comment out if only 1 data set)
path = 'Z:\Lab\Tim\2017-03-10 Racowsky - Human Embs 5 and 2 per O2\2017-03-24 2 per O2\';
file_concat = [path 'ParamsAllSams_samegroupingas5per'  '.xls'];
Tfile = [path 'TimeArrs'];
[num,txt_concat,res_concat]= xlsread(file_concat);
txt_concat(1,:)=[]; res_concat(1,:)=[];
% Concatenate arrays
res = [res; res_concat]; txt = [txt; txt_concat];
dat = load(Tfile);
TimeArrs_concat = dat.TimeArrs; TimeArrs_std_concat = dat.TimeArrs_std;
TimeArrs = [TimeArrs TimeArrs_concat];
TimeArrs_std = [TimeArrs_std TimeArrs_std_concat];
% 
% % Data set 3 (comment out if only 1 data set)
% path = 'Z:\Lab\Tim\2016-10-03 Live Birth Acquisitions\2016-11-10 Batch3\';
% file_concat = [path 'ParamsAllSams'  '.xls'];
% Tfile = [path 'TimeArrs'];
% [num,txt_concat,res_concat]= xlsread(file_concat);
% txt_concat(1,:)=[]; res_concat(1,:)=[];
% res = [res; res_concat]; txt = [txt; txt_concat]; 
% dat = load(Tfile);
% TimeArrs_concat = dat.TimeArrs; TimeArrs_std_concat = dat.TimeArrs_std;
% TimeArrs = [TimeArrs TimeArrs_concat];
% TimeArrs_std = [TimeArrs_std TimeArrs_std_concat];


% Omit rows with a '1' in the 'Exclude' column
txt([res{:,15}]==1,:)=[];
TimeArrs([res{:,15}]==1)=[]; TimeArrs_std([res{:,15}]==1)=[];
res([res{:,15}]==1,:)=[];
% Get the group numbers in an array for easy reference
GrpArr = [res{:,13}];
if isempty(find((~isnan(GrpArr)))) % no groups specified, plot all
    GrpArr = ones(1,length(GrpArr));
    PlotGroups = [1];
    markers = {'o'};
else % User input required
    PlotGroups = [1 2 3];
    %1: All embs
    %2: maybe delayed or dying embs?
    markers = {'o','o','o','o','o','o'};
end

L = size(PlotGroups,2);
r=[1 0 0]; g=[0 1 0]; b=[117 255 255]./255; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1]; p=[160 32 240]./255;
col = [g;b;r;g;b;r];
if AvePlotsBool
    % To better visualize average plots, make individual plots 1/2
    % transparent. Use a 4D color value, 4th value is transparency.
    col = [col 0.5*ones(size(col,1),1)];
end
LegNum = L; % 

AvePlotInd = 1:L;
% AvePlotInd = 1;
%% Calculate average time curves
% Only use time points that are shared by all selected samples
% Find the overlapping time points for this data set so we can calculate an
% average time curve over all samples. Start 'SharedtimePoints' very large
% and take intersect with each sample's 'samtrange' below
for i = AvePlotInd
    SharedTimePoints{i} = 1:10^6;
    Grpinds = find(GrpArr==PlotGroups(i));
    for j = 1:length(Grpinds) % Loop over all samples in this group
        ind = Grpinds(j);
        samtrange = find(TimeArrs(ind).Nbound)&~isnan(TimeArrs(ind).Nbound);
        SharedTimePoints{i} = intersect(find(samtrange),SharedTimePoints{i});
%         if max(TimeArrs(ind).Ntime./60)<85
%             TimeArrs(ind).labels
%         end
    end
end
% Now iterate over time points and average over samples for each time point
for i = AvePlotInd
    Grpinds = find(GrpArr==PlotGroups(i));
    for t = SharedTimePoints{i}
        Ntimes = []; Ftimes = [];
        Nirrs = []; Ntau1s = []; Ntau2s = []; Nbounds = [];
        Firrs = []; Ftau1s = []; Ftau2s = []; Fbounds = [];
        %     redox = Nirr_Z_aves{i}(samtrange)./Firr_Z_aves{i}(samtrange);
        Nirrs_std = []; Ntau1s_std = []; Ntau2s_std = []; Nbounds_std = [];
        Firrs_std = []; Ftau1s_std = []; Ftau2s_std = []; Fbounds_std = [];
        for j = 1:length(Grpinds)
            lab = res(Grpinds(j),1);
            ind = find(strcmp(res(:,1),lab));
            Ntimes = [Ntimes TimeArrs(ind).Ntime(t)];
            Nirrs = [Nirrs TimeArrs(ind).Nirr(t)];
            Nbounds = [Nbounds TimeArrs(ind).Nbound(t)];
            Ntau1s = [Ntau1s TimeArrs(ind).Ntau1(t)];
            Ntau2s = [Ntau2s TimeArrs(ind).Ntau2(t)];
            
            Nirrs_std = [Nirrs_std TimeArrs_std(ind).Nirr(t)];
            Nbounds_std = [Nbounds_std TimeArrs_std(ind).Nbound(t)];
            Ntau1s_std = [Ntau1s_std TimeArrs_std(ind).Ntau1(t)];
            Ntau2s_std = [Ntau2s_std TimeArrs_std(ind).Ntau2(t)];
            
            Ftimes = [Ftimes TimeArrs(ind).Ftime(t)];
            Firrs = [Firrs TimeArrs(ind).Firr(t)];
            Fbounds = [Fbounds TimeArrs(ind).Fbound(t)];
            Ftau1s = [Ftau1s TimeArrs(ind).Ftau1(t)];
            Ftau2s = [Ftau2s TimeArrs(ind).Ftau2(t)];
            
            Firrs_std = [Firrs_std TimeArrs_std(ind).Firr(t)];
            Fbounds_std = [Fbounds_std TimeArrs_std(ind).Fbound(t)];
            Ftau1s_std = [Ftau1s_std TimeArrs_std(ind).Ftau1(t)];
            Ftau2s_std = [Ftau2s_std TimeArrs_std(ind).Ftau2(t)];
        end
        
        AvePlots(i).Ntime(t) = mean(Ntimes); 
        AvePlots(i).Nirr(t) = sum(Nirrs./Nirrs_std.^2)./sum(1./Nirrs_std.^2); 
        AvePlots(i).Ntau1(t) = sum(Ntau1s./Ntau1s_std.^2)./sum(1./Ntau1s_std.^2); 
        AvePlots(i).Ntau2(t) = sum(Ntau2s./Ntau2s_std.^2)./sum(1./Ntau2s_std.^2); 
        AvePlots(i).Nbound(t) = sum(Nbounds./Nbounds_std.^2)./sum(1./Nbounds_std.^2); 
       
        AvePlots(i).Ftime(t) = mean(Ftimes); 
        AvePlots(i).Firr(t) = sum(Firrs./Firrs_std.^2)./sum(1./Firrs_std.^2); 
        AvePlots(i).Ftau1(t) = sum(Ftau1s./Ftau1s_std.^2)./sum(1./Ftau1s_std.^2); 
        AvePlots(i).Ftau2(t) = sum(Ftau2s./Ftau2s_std.^2)./sum(1./Ftau2s_std.^2); 
        AvePlots(i).Fbound(t) = sum(Fbounds./Fbounds_std.^2)./sum(1./Fbounds_std.^2); 
        
        AvePlots(i).redox(t) = AvePlots(i).Nirr(t)/AvePlots(i).Firr(t);
        
        % Note: means are calculated by weighting the sample values by
        % their respective std dev's on the measurement. But the std dev
        % for each param here are calculated by taking simple std dev between samples.
        AvePlots(i).Nirrs_std(t) = std(Nirrs); AvePlots(i).Ntau1s_std(t) = std(Ntau1s); 
        AvePlots(i).Ntau2s_std(t) = std(Ntau2s); AvePlots(i).Nbounds_std(t) = std(Nbounds);
        AvePlots(i).Firrs_std(t) = std(Firrs); AvePlots(i).Ftau1s_std(t) = std(Ftau1s); 
        AvePlots(i).Ftau2s_std(t) = std(Ftau2s); AvePlots(i).Fbounds_std(t) = std(Fbounds);
    end        
end

%% Irrs
close all

figure('position',[200 100 1100 450])

if ~exist('col') col = [(AvePlotInd)'./L (1-(AvePlotInd)'./L) (1-(AvePlotInd)'./L)]; end

specind = [];
for i = AvePlotInd
    Grpinds = find(GrpArr==PlotGroups(i));
    
    subplot(1,2,1)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ntime./60; 
        y = TimeArrs(ind).Nirr;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('NADH Brightness')
    
    subplot(1,2,2)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ftime./60;
        y = TimeArrs(ind).Firr;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('FAD Brightness')
        
end
if AvePlotsBool
    for i = AvePlotInd
        non0ind = find(AvePlots(i).Nirr);
        subplot(1,2,1)
        h(i) = plot(AvePlots(i).Ntime(non0ind)./60,AvePlots(i).Nirr(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
        subplot(1,2,2)
        h(i) = plot(AvePlots(i).Ftime(non0ind)./60,AvePlots(i).Firr(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
    end
end
if exist('GroupNames') legend(h,GroupNames(1:LegNum),'location','best');end 
clear h;
set(gcf,'PaperPositionMode','auto')
set(gcf,'color',[1 1 1])
saveas(gcf,[OutFile '_irrs.jpg'],'jpg')
saveas(gcf,[OutFile '_irrs.fig'],'fig')

%% Fbound
close all

figure('position',[200 100 1100 450])

specind = [];
for i = AvePlotInd
    Grpinds = find(GrpArr==PlotGroups(i));
    subplot(1,2,1)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ntime./60;
        y = TimeArrs(ind).Nbound;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('NADH Fraction Bound')
    
    subplot(1,2,2)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ftime./60;
        y = TimeArrs(ind).Fbound;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('FAD Fraction Bound')
end
if AvePlotsBool
    for i = AvePlotInd
        non0ind = find(AvePlots(i).Nirr);
        subplot(1,2,1)
        h(i) = plot(AvePlots(i).Ntime(non0ind)./60,AvePlots(i).Nbound(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
        subplot(1,2,2)
        h(i) = plot(AvePlots(i).Ftime(non0ind)./60,AvePlots(i).Fbound(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
    end
end
if exist('GroupNames') legend(h,GroupNames(1:LegNum),'location','best');end 
set(gcf,'PaperPositionMode','auto')
set(gcf,'color',[1 1 1])
saveas(gcf,[OutFile '_bounds.jpg'],'jpg')
saveas(gcf,[OutFile '_bounds.fig'],'fig')

%% Tau1
close all

figure('position',[200 100 1100 450])


specind = [];
for i = AvePlotInd
    Grpinds = find(GrpArr==PlotGroups(i));
    subplot(1,2,1)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ntime./60;
        y = TimeArrs(ind).Ntau1;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('NADH \tau_{1}')
    
    subplot(1,2,2)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ftime./60;
        y = TimeArrs(ind).Ftau1;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('FAD \tau_{1}')
end
if AvePlotsBool
    for i = AvePlotInd
        non0ind = find(AvePlots(i).Nirr);
        subplot(1,2,1)
        h(i) = plot(AvePlots(i).Ntime(non0ind)./60,AvePlots(i).Ntau1(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
        subplot(1,2,2)
        h(i) = plot(AvePlots(i).Ftime(non0ind)./60,AvePlots(i).Ftau1(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
    end
end
if exist('GroupNames') legend(h,GroupNames(1:LegNum),'location','best');end 
clear h
set(gcf,'PaperPositionMode','auto')
set(gcf,'color',[1 1 1])
saveas(gcf,[OutFile '_tau1s.jpg'],'jpg')
saveas(gcf,[OutFile '_tau1s.fig'],'fig')

%% Tau2
close all

figure('position',[200 100 1100 450])


specind = [];
for i = AvePlotInd
    Grpinds = find(GrpArr==PlotGroups(i));
    subplot(1,2,1)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ntime./60;
        y = TimeArrs(ind).Ntau2;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('NADH \tau_{2}')
    
    subplot(1,2,2)
    for j = 1:length(Grpinds)
        ind = Grpinds(j);
        x = TimeArrs(ind).Ftime./60;
        y = TimeArrs(ind).Ftau2;
        non0indtmp = find(TimeArrs(ind).Nirr~=0&~isnan(TimeArrs(ind).Nirr));
        x = x(non0indtmp); y = y(non0indtmp);
        plot(x,y,'color',col(PlotGroups(i),:)); hold on
    end
    xlabel('time (h)'); ylabel('FAD \tau_{2}')
end
if AvePlotsBool
    for i = AvePlotInd
        non0ind = find(AvePlots(i).Nirr);
        subplot(1,2,1)
        h(i) = plot(AvePlots(i).Ntime(non0ind)./60,AvePlots(i).Ntau2(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
        subplot(1,2,2)
        h(i) = plot(AvePlots(i).Ftime(non0ind)./60,AvePlots(i).Ftau2(non0ind),'color',col(PlotGroups(i),1:3),'linewidth',3);
    end
end
if exist('GroupNames') legend(h,GroupNames(1:LegNum),'location','best');end
set(gcf,'PaperPositionMode','auto')
set(gcf,'color',[1 1 1])
saveas(gcf,[OutFile '_tau2s.jpg'],'jpg')
saveas(gcf,[OutFile '_tau2s.fig'],'fig')