function FLIMBatchParamPlot_wROIMovie(names,wrtpath,MaskLab)
% Plots params from a saved data set from FLIM analysis (FittingGUI_ver3_4)

% clear all;
% % names = {'Z:\Lab\Tim\2016-10-03 Live Birth Acquisitions\2017-02-21 Batch4\s1_a1\fits_Pos0_mask1_fxshft.mat'};
% names = {'C:\Users\Tim\Documents\Academic - Research\Data\2017-03-10 Racowsky\s1_a1_HumEmbs\fits_Pos0_BlastOnly_fxshft.mat'};
% wrtpath = 'C:\Users\Tim\Documents\Academic - Research\Data\2017-03-10 Racowsky\s1_a1_HumEmbs\PlotswMovies\fits_Pos0_BlastOnly_fxshft';
% MaskLab = 'BlastOnly';

if ~exist('MaskLab') MaskLab = 'Masks'; end

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
% 2017-04-07: Changed graph styles to be current with recent batch plot
% 2015-10-15: Previous version stepped through nameinds. Change to walk
% through time points and display a merge of all z slices.
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

[a,b] = mkdir(wrtpath);

%% Group data in one cell. Also construct a cell that can be written in spreadsheet form
AllDecays = [];
AllFNames = {};
% Spr{3,1} = 'shift, shift_par';
% Spr{4,1} = 'shift_perp';
% Spr{5,1} = 'bg_par';
% Spr{6,1} = 'bg_perp';
% Spr{7,1} = 'tau1';
% Spr{8,1} = 'r0_1';
% Spr{9,1} = 'theta1';
% Spr{10,1} = 'f';
% Spr{11,1} = 'tau2';
% Spr{12,1} = 'r0_2';
% Spr{13,1} = 'theta2';
% Spr{14,1} = 'Chi2_red';

tm = [];
Chs = [];
Ims = {};
imw = 0; imh = 0;
Allmasks = {};
NmTinds = [];
tstructs = [];
for i=1:length(names)
    % Times and channel info
    load(names{i});
    load([UpOneDir(names{i}) 'name_indexes'])
    decbools = ~cellfun('isempty',decays_fits_struct);
    % Load images and masks, draw ROIs and crop, and store in a cell to
    % later make a movie next to the params
    slashes = strfind(names{i},'\');
    fname = names{i}(slashes(end)+1:end);
    dashes =  strfind(fname,'_');
    run = names{i}(slashes(end-1)+1:slashes(end)-1);
    posind = strfind(fname,'Pos');
    NextDash = find(dashes>posind); NextDash = dashes(NextDash(1));
    uManPos = fname(posind+3:NextDash-1);
    
    % Deterimine whether there are individual embyro masks, or only one
    % single mask for the position. If so, there will be a label, 'mask' in
    % the fit file name, followed by the mask #
    maskind = strfind(fname,'mask');
    if ~isempty(maskind)
        NextDash = find(dashes>maskind); NextDash = dashes(NextDash(1));
        masknum = num2str(fname(maskind+4:NextDash-1));
    else
        masknum = 1; % Only one single mask.
    end
    %
    %     strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
    %     uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    %     maskind = strfind(names{i},'mask');
    %     masknum = names{i}(maskind+4:end);
    %     strnums = sscanf(masknum,'%g'); %Find the numbers in the name
    %     if isempty(strnums)
    %         load([UpOneDir(names{i}) run '_Pos' num2str(uManPos) '_GenMask.mat'])
    %         masknum = 1;
    %     else
    %         masknum = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    %         load([UpOneDir(names{i}) run '_Pos' num2str(uManPos) '_Masks.mat'])
    %     end
    
    load([UpOneDir(names{i}) run '_Pos' uManPos '_' MaskLab '.mat'])
    
    % Fill cells, one time point at at time
    trg = unique(cell2mat(nameinds([nameinds{:,7}]>-1,3)));
    clear tstruct
    for j = 1:length(trg)
        currTbools = strcmp(nameinds(:,2),uManPos)&([nameinds{:,3}]==trg(j))';%&([nameinds{:,7}]>-1))';
        currTinds = find(currTbools);
        if isempty(decays_fits_struct(currTbools&decbools))
            continue;
        end
        % Try to use the middle image in the z-stack if it was found. If
        % not, try for the highest, then the lowest. This is a Wonky way,
        % but I couldn't figure out something more elegant. It works.
        clear ind
        for k = 1:length(currTinds)
            if [nameinds{currTinds(k),5}]==1&[nameinds{currTinds(k),7}]>-1&~isempty(Masks{nameinds{currTinds(k),6}}(masknum).NL)
                ind = k;
                break;
            end
        end
        if ~exist('ind')
            for k = 1:length(currTinds)
                if [nameinds{currTinds(k),5}]==2&[nameinds{currTinds(k),7}]>-1&~isempty(Masks{nameinds{currTinds(k),6}}(masknum).NL)
                    ind = k;
                    break;
                end
            end
        end
        if ~exist('ind')
            for k = 1:length(currTinds)
                if [nameinds{currTinds(k),5}]==0&[nameinds{currTinds(k),7}]>-1&~isempty(Masks{nameinds{currTinds(k),6}}(masknum).NL)
                    ind = k;
                    break;
                end
            end
        end
        if ~exist('ind') error('No mask found. Something went wrong.'); end
        
        tstruct(j).mask = Masks{nameinds{currTinds(ind),6}}(masknum).NLper;
        tstruct(j).im = imread([UpOneDir(names{i}) 'sorted_sdts\IntTiffs_Pos' num2str(uManPos) '_' run '\fr' num2str(nameinds{currTinds(ind),6},'%05i') '.tif']);
        
        
        % Correspondence index between struct and AllDecays
        l = length([nameinds{currTbools&decbools,7}]');
        NmTind = [ones(l,1).*i ones(l,1).*j (1:l)'];
        
        [y,x] = find(tstruct(j).mask);
        if (max(x)-min(x))>imw imw = (max(x)-min(x)); end
        if (max(y)-min(y))>imh imh = (max(y)-min(y)); end
        AllDecays = [AllDecays; decays_fits_struct(currTbools&decbools)];
        tm = [tm; [nameinds{currTinds(ind),7}]'; [nameinds{currTinds(ind+floor(size(currTinds,1)/2)),7}]'];
        Chs = [Chs; nameinds(currTinds(ind),4); nameinds(currTinds(ind+floor(size(currTinds,1)/2)),4)];
        %         tm = [tm; [nameinds{currTbools&decbools,7}]'];
        %         Chs = [Chs; nameinds(currTbools&decbools,4)];
        NmTinds = [NmTinds; NmTind];
    end
    tstructs = [tstructs tstruct];
end
xdim = size(tstruct(end).mask,2); ydim = size(tstruct(end).mask,1);

tm = (tm - min(tm))*86400/timeunit;

L = length(AllDecays);
for i = 1:L
    fitres = AllDecays{i}.fit_result;
    Chisqred(i) = AllDecays{i}.Chi_sq;
    AllParams{i} = fitres;
    irrs(i) = AllDecays{i}.irr;
    numphot(i) = sum(AllDecays{i}.decay);
    %     Spr{2,2+(i-1)} = AllDecays{i}.name;
    %     Spr{3,2+(i-1)} = fitres(1,1);
    %     Spr{4,2+(i-1)} = fitres(2,1);
    %     Spr{5,2+(i-1)} = fitres(3,1);
    %     Spr{6,2+(i-1)} = fitres(4,1);
    %     Spr{7,2+(i-1)} = fitres(5,1);
    %     Spr{8,2+(i-1)} = fitres(6,1);
    %     Spr{9,2+(i-1)} = fitres(7,1);
    %     Spr{10,2+(i-1)} = fitres(8,1);
    %     Spr{11,2+(i-1)} = fitres(9,1);
    %     Spr{12,2+(i-1)} = fitres(10,1);
    %     Spr{13,2+(i-1)} = fitres(11,1);
    %     Spr{14,2+(i-1)} = Chisqred(i);
end
ParamsArr = reshape(cell2mat(AllParams),5,3,length(AllDecays));
% Define plot arrays
for i = 1:size(ParamsArr,3)
    f(i) = reshape(ParamsArr(4,1,i),1,1);
    shifts(i) = reshape(ParamsArr(1,1,i),1,1);
    tau1(i) = reshape(ParamsArr(3,1,i),1,1);
    tau2(i) = tau1(i)*reshape(ParamsArr(5,1,i),1,1);
end

Nind = strcmp(Chs,'NADH');
Find = strcmp(Chs,'FAD');

%% Put into spreadsheet
% xlswrite([pathname [names{:}] '.xls'],Spr)

%% If no arguments entered, plot in one big plot.
% if nargin<2
fh = gcf;%% If no arguments entered, plot in one big plot.
% if nargin<2
fh = gcf;
clf(fh)
set(fh,'Name','Parameters Comparison','Units','normalized','Position',[.1 .1 .8 .8],'Color','w'); hold on;

% Identify bad fits
misfits = [];
misfits = f==0|f==1; misfits = misfits';
if isempty(find(Nind&~misfits)) misfits(Nind) = 0; end % plot error below if no points present
if isempty(find(Find&~misfits)) misfits(Find) = 0; end


% SHIFTS and CHI^2's
splts(1) = subplot(2,3,1,'position',[0.07 0.55 0.27 0.38]); hold on;
title('Shifts (o) and \chi^2s (.)','fontsize',12)
grid on;
if ~isempty(find(Nind))&isempty(find(Find))
    [ax_chi,h1,h2] =plotyy(tm(Nind),shifts(Nind),tm(Nind),Chisqred(Nind));
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
elseif ~isempty(find(Find))&isempty(find(Nind))
    [ax_chi,h1,h2] =plotyy(tm(Find),shifts(Find),tm(Find),Chisqred(Find));
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    xlim([min(tm) max(tm)])
elseif ~isempty(find(Nind))&~isempty(find(Find))
    [ax_chi,h1,h2] =plotyy(tm(Nind),shifts(Nind),tm(Nind),Chisqred(Nind));
    ax_chi(1).XLim = [min(tm) max(tm)];
    ax_chi(2).XLim = [min(tm) max(tm)];
    ax_chi(1).YColor = [0 0 0];
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    hold(ax_chi(1),'on');
    plot(ax_chi(1),tm(Find),shifts(Find),'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    ax_chi(1).YLim(1) =  min(shifts); ax_chi(1).YLim(2) =  max(shifts);
    hold(ax_chi(2),'on');
    plot(ax_chi(2),tm(Find),Chisqred(Find),'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    ylim(ax_chi(2),[min(Chisqred)*.95 max(Chisqred)*1.05])
    ax_chi(1).FontSize = 11; ax_chi(2).FontSize = 11;
    ax_chi(1).LabelFontSizeMultiplier=1.2; ax_chi(2).LabelFontSizeMultiplier=1.2;
    % Set reasonable ticks
    ndig = ceil(log10(abs(max(Chisqred)-min(Chisqred))));
    ytk1 = round(min(shifts./10^(ndig-1))); ytk1 = ytk1*10^(ndig-1); ytk1 = ytk1(1); ytk2 = round(max(shifts./10^(ndig-1))); ytk2 = ytk2*10^(ndig-1); ytk2 = ytk2(1);
    set(ax_chi(1),'YTick',ytk1:(ytk2-ytk1)/5:ytk2)
    xlim([min(tm) max(tm)])
    % Set reasonable ticks
    ndig = ceil(log10(abs(max(Chisqred)-min(Chisqred))));
    ytk1 = round(min(Chisqred./10^(ndig-1))); ytk1 = ytk1*10^(ndig-1); ytk1 = ytk1(1); ytk2 = round(max(Chisqred./10^(ndig-1))); ytk2 = ytk2*10^(ndig-1); ytk2 = ytk2(1);
    set(ax_chi(2),'YTick',ytk1:(ytk2-ytk1)/5:ytk2)
end
ax_chi(1).YLabel.String =  'Shifts';
% ax_chi(2).YLabel.String =  '\chi^2s';
grid on;
set(gca,'fontsize',12)


% FBOUNDS
splts(2) = subplot(2,3,2,'position',[0.4 0.55 0.27 0.38]);
axes(splts(2)); hold on;
title('Population Fractions Bound','fontsize',12)
grid on; hold on;
if ~isempty(find(Nind))&isempty(find(Find))
    plot(tm(Nind&~misfits),f(Nind&~misfits),'ob','MarkerSize',11,'LineWidth',2);
elseif ~isempty(find(Find))&isempty(find(Nind))
    plot(tm(Find&~misfits),f(Find&~misfits),'og','MarkerSize',11,'LineWidth',2);
elseif ~isempty(find(Nind))&~isempty(find(Find))
    [ax_f,h1,h2] =plotyy(tm(Nind&~misfits),f(Nind&~misfits),tm(Find&~misfits),1-f(Find&~misfits));
    ax_f(1).YColor = [0 0 1]; ax_f(2).YColor = [0 1 0];
    axes(ax_f(2)); hold on;
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
%     ax_f(1).YLabel.String =  'NADH fbound';
%     ax_f(2).YLabel.String =  'FAD fbound';
    ax_f(1).FontSize = 11; ax_f(2).FontSize = 11;
    ax_f(1).LabelFontSizeMultiplier=1.2; ax_f(2).LabelFontSizeMultiplier=1.2;
    %     Set reasonable ticks and limits
    GoodPlotyyTicks(ax_f,tm,f(Nind&~misfits),1-f(Find&~misfits));
end


% LIFETIMES
splts(4) = subplot(2,3,4,'position',[0.07 0.07 0.27 0.38]); hold on;
title('Lifetimes (ns): short (o) and long (\Delta)','fontsize',12)
grid on;
twoexpB = 0;
if ~isempty(find(Nind))&isempty(find(Find))
    [ax_life,h1,h2] =plotyy(tm(Nind&~misfits),tau1(Nind&~misfits),tm(Nind&~misfits),tau2(Nind&~misfits));
    set(h1,'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
elseif ~isempty(find(Find))&isempty(find(Nind))
    [ax_life2,h1,h2] =plotyy(tm(Find&~misfits),tau1(Find&~misfits),tm(Find&~misfits),tau2(Find&~misfits));
    set(h1,'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    xlim([min(tm) max(tm)])
elseif ~isempty(find(Nind))&~isempty(find(Find))
    % Tau1 is lifetime, tau2 is short
    [ax_life,h1,h2] =plotyy(tm(Nind&~misfits),tau2(Nind&~misfits),tm(Nind&~misfits),tau1(Nind&~misfits));
    ax_life(1).XLim = [min(tm) max(tm)];
    ax_life(2).XLim = [min(tm) max(tm)];
    ax_life(1).YColor = [0 0 0];
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    hold(ax_life(1),'on');
    plot(ax_life(1),tm(Find&~misfits),tau2(Find&~misfits),'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    ax_life(1).YLim = [min(tau2(~misfits))*.99 max(tau2(~misfits))*1.01];
    hold(ax_life(2),'on');
    plot(ax_life(2),tm(Find),tau1(Find&~misfits),'LineStyle','none','Marker','^','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    ax_life(2).YLim = [min(tau1(~misfits))*.99 max(tau1(~misfits))*1.01];
    ax_life(1).FontSize = 11; ax_life(2).FontSize = 11;
    ax_life(1).LabelFontSizeMultiplier=1.2; ax_life(2).LabelFontSizeMultiplier=1.2;
    %     Set reasonable ticks and limits
    GoodPlotyyTicks(ax_life,tm,tau2(~misfits),tau1(~misfits));
end
% ax_life(1).YLabel.String =  'Short Lifetimes (o)';
% ax_life(2).YLabel.String =  'Long Lifetimes (\Delta)';
xlabel(tlab,'fontsize',12)


% IRRADIANCES
splts(5) = subplot(2,3,5,'position',[0.4 0.07 0.27 0.38]); hold on;
title('Irradiances (phots/pixel/scan)','fontsize',12)
grid on;
if ~isempty(find(Nind))&isempty(find(Find))
    [ax_irr,h1,h2] =plotyy(tm(Nind),irrs(Nind),tm(Nind),numphot(Nind));
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
elseif ~isempty(find(Find))&isempty(find(Nind))
    [ax_irr2,h1,h2] =plotyy(tm(Find),irrs(Find),tm(Find),numphot(Find));
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    set(h2,'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    xlim([min(tm) max(tm)])
elseif ~isempty(find(Nind))&~isempty(find(Find))
    [ax_irr,h1,h2] =plotyy(tm(Nind),irrs(Nind),tm(Find),irrs(Find));
    ax_irr(1).XLim = [min(tm) max(tm)];
    ax_irr(2).XLim = [min(tm) max(tm)];
    ax_irr(1).YColor = [0 0 1]; ax_irr(2).YColor = [0 1 0];
    set(h1,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','b','LineWidth',2);
    set(h2,'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','g','LineWidth',2);
    hold(ax_irr(1),'on');
    plot(ax_irr(1),tm(Nind&misfits),irrs(Nind&misfits),'LineStyle','none','Marker','x','Markersize',11,'MarkerEdgeColor','r','LineWidth',2);
    hold(ax_irr(2),'on');
    plot(ax_irr(2),tm(Find&misfits),irrs(Find&misfits),'LineStyle','none','Marker','x','Markersize',11,'MarkerEdgeColor','r','LineWidth',2);
    ax_irr(1).YLim = [min(irrs(Nind))*0.99 max(irrs(Nind))*1.01];
    ax_irr(2).YLim = [min(irrs(Find))*0.99 max(irrs(Find))*1.01];
    ax_irr(1).FontSize = 11; ax_irr(2).FontSize = 11;
    ax_irr(1).LabelFontSizeMultiplier=1.2; ax_irr(2).LabelFontSizeMultiplier=1.2;
    %     Set reasonable ticks and limits
    GoodPlotyyTicks(ax_irr,tm,irrs(Nind),irrs(Find));
end
grid on;
% ax_irr(1).YLabel.String =  'NADH Irr';
% ax_irr(2).YLabel.String =  'FAD Irr';
xlabel(tlab,'fontsize',12)
% legend(h1,'Irr','location','nw')
% legend(h2,'Photons','location','ne')

set(gcf,'PaperPositionMode','auto')
%% Show Ims and make a movie
% Do every other frame. Since only frames with NADH and FAD are accepted,
% the channels are exactly alternating at this point. Therefore do every
% other frame, so you can highlight both the NADH and FAD single points at
% once
border = 10;

% Adjust NmTinds to have linear time instead of starting over with each run
NmTinds2=NmTinds;
for i = 1:(length(names)-1)
    NmTinds2(NmTinds(:,1)==i+1,2) = NmTinds2(NmTinds(:,1)==i+1,2) + max(NmTinds2(NmTinds(:,1)==i,2));
end

for i = 1:2:size(tm,1)
    mask = [zeros(ydim+2*imh,imw) [zeros(imh,xdim);tstructs(NmTinds2(i,2)).mask;zeros(imh,xdim)] zeros(ydim+2*imh,imw)];
    Nim = [zeros(ydim+2*imh,imw) [zeros(imh,xdim);tstructs(NmTinds2(i,2)).im(:,1:xdim);zeros(imh,xdim)] zeros(ydim+2*imh,imw)];
    Fim = [zeros(ydim+2*imh,imw) [zeros(imh,xdim);tstructs(NmTinds2(i,2)).im(:,xdim+1:2*xdim);zeros(imh,xdim)] zeros(ydim+2*imh,imw)];
    [y,x] = find(mask); CoM = [mean(x) mean(y)];
    % CROP
    CropNim = Nim(round(CoM(2)-imh/2-border):round(CoM(2)+imh/2+border),round(CoM(1)-imw/2-border):round(CoM(1)+imw/2+border));
    CropFim = Fim(round(CoM(2)-imh/2-border):round(CoM(2)+imh/2+border),round(CoM(1)-imw/2-border):round(CoM(1)+imw/2+border));
    mask = mask(round(CoM(2)-imh/2-border):round(CoM(2)+imh/2+border),round(CoM(1)-imw/2-border):round(CoM(1)+imw/2+border));
    [y,x] = find(mask); CropCoords = [x,y];
    IntInds = find(mask);
    
    % NADH
    sp(1) = plot(ax_chi(1),tm(i),shifts(i),'oc','MarkerSize',11,'LineWidth',2);
    sp(5) = plot(ax_chi(2),tm(i),Chisqred(i),'.c','MarkerSize',11,'LineWidth',2);
    sp(2) = plot(ax_f(1),tm(i),f(i),'oc','MarkerSize',11,'LineWidth',2);
    sp(3) = plot(ax_life(2),tm(i),tau1(i),'^c','MarkerSize',11,'LineWidth',2);
    sp(4) = plot(ax_life(1),tm(i),tau2(i),'oc','MarkerSize',11,'LineWidth',2);
    sp(6) = plot(ax_irr(1),tm(i),irrs(i),'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','c','LineWidth',2);
    %     sp(7) = plot(ax(2),tm(i),numphot(i),'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','c','LineWidth',2);
    
    % FAD
    sp(8) = plot(ax_chi(1),tm(i+1),shifts(i+1),'om','MarkerSize',11,'LineWidth',2);
    sp(12) = plot(ax_chi(2),tm(i+1),Chisqred(i+1),'.m','MarkerSize',11,'LineWidth',2);
    sp(9) = plot(ax_f(2),tm(i+1),1-f(i+1),'om','MarkerSize',11,'LineWidth',2);
    sp(10) = plot(ax_life(2),tm(i+1),tau1(i+1),'^m','MarkerSize',11,'LineWidth',2);
    sp(11) = plot(ax_life(1),tm(i+1),tau2(i+1),'om','MarkerSize',11,'LineWidth',2);
    sp(13) = plot(ax_irr(2),tm(i+1),irrs(i+1),'LineStyle','none','Marker','o','Markersize',11,'MarkerEdgeColor','m','LineWidth',2);
    %     sp(14) = plot(ax(2),tm(i+1),numphot(i+1),'LineStyle','none','Marker','.','Markersize',11,'MarkerEdgeColor','m','LineWidth',2);
    
    % Imgs
    subplot(2,3,3)
    Ints = CropNim(IntInds);
    Nih = imshow(CropNim,[min(Ints) max(Ints)]); hold on;
    set(gca,'position',[0.7 0.55 0.27 0.38])
    Nc = plot(CropCoords(:,1),CropCoords(:,2),'.b'); hold off;
    subplot(2,3,6)
    Ints = CropFim(IntInds);
    Fih = imshow(CropFim,[min(Ints) max(Ints)]); hold on;
    set(gca,'position',[0.7 0.07 0.27 0.38])
    Fc = plot(CropCoords(:,1),CropCoords(:,2),'.g'); hold off;
    
    % Save, then delete the single point plots
    set(splts(1),'position',[0.07 0.55 0.27 0.38]); %compensation for total matlab glitch
    set(splts(2),'position',[0.4 0.55 0.27 0.38]); 
    set(splts(4),'position',[0.07 0.07 0.27 0.38]);
    set(splts(5),'position',[0.4 0.07 0.27 0.38]);
    
    saveas(gcf,[wrtpath '\T' num2str(tm(i),3) tlab(end-2:end) '.jpg'],'jpg')
    delete(sp);
end


