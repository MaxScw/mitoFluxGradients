
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

% path = pwd;
% path = UpOneDir(path);

close all
% figure('units','normalized','position',[.1 .1 .8 .8])
figure('position',[200 100 1300 850])

%% USER INPUT
OutFile = ['Z:\Lab\Tim\2017-03-10 Racowsky - Human Embs 5 and 2 per O2\Analysis and Plots\AllEmbsDiffCol_alltime'];

% If text labels for all data points are desired (they are entered on the
% spreadsheet)
TextBool = 0;

% Data set 1
path = 'Z:\Lab\Tim\2017-03-10 Racowsky - Human Embs 5 and 2 per O2\';
file = [path 'ParamsAllSams_alltime'  '.xls'];
file_std = [path 'ParamsAllSams_alltime'  '_std.xls'];
[num,txt,res]= xlsread(file);
[num_std,txt_std,res_std]= xlsread(file_std);
% Remove 1st row of headers
txt(1,:)=[]; res(1,:)=[]; res_std(1,:)=[];

% Data set 2 (comment out if only 1 data set)
path = 'Z:\Lab\Tim\2017-03-10 Racowsky - Human Embs 5 and 2 per O2\2017-03-24 2 per O2\';
file_concat = [path 'ParamsAllSams_alltime'  '.xls'];
[num,txt_concat,res_concat]= xlsread(file_concat);
file_std_concat = [path 'ParamsAllSams_alltime'  '_std.xls'];
[num,txt_std_concat,res_std_concat]= xlsread(file_std_concat);
txt_concat(1,:)=[]; res_concat(1,:)=[]; res_std_concat(1,:)=[];
% Concatenate arrays
res = [res; res_concat]; txt = [txt; txt_concat];
res_std = [res_std; res_std_concat];
% 
% % Data set 2 (comment out if only 1 data set)
% path = 'Z:\Lab\Tim\2016-10-03 Live Birth Acquisitions\2016-11-10 Batch3\';
% file_concat = [path 'ParamsAllSams'  '.xls'];
% [num,txt_concat,res_concat]= xlsread(file_concat);
% file_std_concat = [path 'ParamsAllSams_std'  '.xls'];
% [num,txt_std_concat,res_std_concat]= xlsread(file_std_concat);
% txt_concat(1,:)=[]; res_concat(1,:)=[]; res_std_concat(1,:)=[];
% res = [res; res_concat]; txt = [txt; txt_concat];
% res_std = [res_std; res_std_concat];


% Omit rows with a '1' in the 'Exclude' column
txt([res{:,15}]==1,:)=[];
res([res{:,15}]==1,:)=[];
% Get the group numbers in an array for easy reference
GrpArr = [res{:,13}];
if isempty(find((~isnan(GrpArr)))) % no groups specified, plot all
    GrpArr = ones(1,length(GrpArr));
    PlotGroups = [1];
    markers = {'o'};
else % User input required
    PlotGroups = [1 2 3 4 5 6];
    %1: All embs
    %2: maybe delayed or dying embs?
    markers = {'o','o','o','x','x','x'};
end

L = size(PlotGroups,2);
r=[1 0 0]; g=[0 1 0]; b=[117 255 255]./255; o = [.9 0.45 0]; k = [0 0 0]; w = [1 1 1]; p=[160 32 240]./255;
col = [g;b;r;g;b;r];
% col = [(1:L)'./L (1-(1:L)'./L) (1-(1:L)'./L)];
LegNum = L; % 



%% PLOTs

for i = 1:L
    PlotInds = find(GrpArr==PlotGroups(i));
    if ~isnan(res{PlotInds(1),14}) GroupNames{i} = res{PlotInds(1),14}; end
    
    
    % Fracs bound
    subplot(2,2,1)
    h = plot(cell2mat(res(PlotInds,7)),cell2mat(res(PlotInds,10)),markers{PlotGroups(i)},'markersize',10,'markerfacecolor',col(PlotGroups(i),:),'color',col(PlotGroups(i),:));hold on
    
    % Tau 1's
    subplot(2,2,2)
    h = plot(cell2mat(res(PlotInds,5)),cell2mat(res(PlotInds,8)),markers{PlotGroups(i)},'markersize',10,'markerfacecolor',col(PlotGroups(i),:),'color',col(PlotGroups(i),:));hold on
    
    % Tau 2's
    subplot(2,2,3)
    h = plot(cell2mat(res(PlotInds,6)),cell2mat(res(PlotInds,9)),markers{PlotGroups(i)},'markersize',10,'markerfacecolor',col(PlotGroups(i),:),'color',col(PlotGroups(i),:));hold on

    % Irrs
    subplot(2,2,4)
    h = plot(cell2mat(res(PlotInds,2)),cell2mat(res(PlotInds,3)),markers{PlotGroups(i)},'markersize',10,'markerfacecolor',col(PlotGroups(i),:),'color',col(PlotGroups(i),:));hold on
end

if TextBool
    % % For automatically getting text labels from names
    % for i = [specrg{:}]
    %     tind1 = strfind(res{i,1},'em')+2; tind2 = strfind(res{i,1},'P')-1;
    %     texts{i} = res{i,1}(tind1:tind2);
    % end
    texts = res(:,16);
    for i = 1:L
        PlotInds = find(GrpArr==PlotGroups(i)&~strcmp(txt(:,16),'')');
        
        % Fracs bound
        subplot(2,2,1)
        text(cell2mat(res(PlotInds,7)),cell2mat(res(PlotInds,10)),texts(PlotInds))
        
        % Tau 1's
        subplot(2,2,2)
        text(cell2mat(res(PlotInds,5)),cell2mat(res(PlotInds,8)),texts(PlotInds))
        
        % Tau 2's
        subplot(2,2,3)
        text(cell2mat(res(PlotInds,6)),cell2mat(res(PlotInds,9)),texts(PlotInds))
        
        % Irrs
        subplot(2,2,4)
        text(cell2mat(res(PlotInds,2)),cell2mat(res(PlotInds,3)),texts(PlotInds))
    end
end


% Formatting and legends
subplot(2,2,1)
% legend(GroupNames(1:LegNum),'location','best')
xlabel('NADH frac bound','fontsize',12); ylabel('FAD frac bound','fontsize',12)
xl = get(gca,'xLim'); yl = get(gca,'yLim');
w = abs(xl(1)-xl(2));xlim([xl(1)-w./10 xl(2)+w./10]);
h = abs(yl(1)-yl(2));ylim([yl(1)-h./10 yl(2)+h./10]);

subplot(2,2,2)
% legend(GroupNames(1:LegNum),'location','best')
xlabel('NADH \tau_1','fontsize',12); ylabel('FAD \tau_1','fontsize',12)
xl = get(gca,'xLim'); yl = get(gca,'yLim');
w = abs(xl(1)-xl(2));xlim([xl(1)-w./10 xl(2)+w./10]);
h = abs(yl(1)-yl(2));ylim([yl(1)-h./10 yl(2)+h./10]);

subplot(2,2,3)
% legend(GroupNames(1:LegNum),'location','best')
xlabel('NADH \tau_2','fontsize',12); ylabel('FAD \tau_2','fontsize',12)
xl = get(gca,'xLim'); yl = get(gca,'yLim');
set(gca,'fontsize',11)
w = abs(xl(1)-xl(2));xlim([xl(1)-w./10 xl(2)+w./10]);
h = abs(yl(1)-yl(2));ylim([yl(1)-h./10 yl(2)+h./10]);

subplot(2,2,4)
xlabel('NADH Irr','fontsize',12); ylabel('FAD Irr','fontsize',12)
xl = get(gca,'xLim'); yl = get(gca,'yLim');
set(gca,'fontsize',11)
w = abs(xl(1)-xl(2));xlim([xl(1)-w./10 xl(2)+w./10]);
h = abs(yl(1)-yl(2));ylim([yl(1)-h./10 yl(2)+h./10]);
if exist('GroupNames') legend(GroupNames(1:LegNum),'location','best');end


% Save
set(gcf,'PaperPositionMode','auto')
% hgexport(gcf,[OutFile '.jpg'])
set(gcf,'color',[1 1 1])
saveas(gcf,[OutFile '.jpg'],'jpg')
saveas(gcf,[OutFile '.fig'],'fig')
% d = get(gcf,'Position');
% screencapture(gcf,[d(3)*.07, d(4)*.05, d(3)*.85, d(4)*.9],[OutFile '.jpg']);
