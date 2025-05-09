function GoodPlotyyTicks(axs,minticks)%(axs,ydata1,ydata2,minticks)
% Function for efficiently setting tight limits and TickLabels for plots
% created with plotyy.
% inputs:
% -axs: axsis handle of the plotyy plot, ie, '[axs,h1,h2] =plotyy(...'
% -xdata: x-axsis data, assumed to be the same for both ydata sets.
% -ydata1and2: y-axsis data array sets used in the plotting
% -minticks: optional number of ticks

% To use on a 1-y plot, just enter the handle for a 1-y axsis and ydata1, only

% clear all; close all;
% uiopen('C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\test.fig',1)
% axs = gca;
% ydata1 = get(gca,'YLim');

% Account for cases where one of the data sets is absent.
% if length(ydata1)<2 ydata1=ydata2; end
% if exist('ydata2')&ydata2>-1 if length(ydata2)<2 ydata2=ydata1; end; end
if ~exist('axs') axs=gca; end
if ~exist('minticks') minticks=4; end


for anum = 1:length(axs)
    Tcks = get(axs(anum),'YTick');
    YL = get(axs(anum),'YLim');
    if length(Tcks)<minticks
        L = [min(Tcks) max(Tcks)];
        Lgap = YL(2)-YL(1); % gap calculation
        bf = .01;
        while length(Tcks)<minticks & bf<.5
            set(gca,'YLim',[YL(1)-Lgap*bf YL(2)+Lgap*bf]);
            Tcks = get(axs(anum),'YTick');
            bf = bf + 0.01;
        end
        YL = get(axs(anum),'YLim');
        Lgap = YL(2)-YL(1); 
        set(axs(anum),'YLim',[mean(Tcks)-Lgap/1.9 mean(Tcks)+Lgap/1.9])
    end
end



% % Old method... buggy
% FirstNonEqualSigFig = log10(10.^floor(log10(abs(Lgap))));
% if FirstNonEqualSigFig<-2 FirstNonEqualSigFig = -2; end
% if Lgap > 0
%     % Recalculate cap to be to nearest sig fig after leading
%     L = round(L*10^-(FirstNonEqualSigFig-1))*10^(FirstNonEqualSigFig-1);
%     Lgap = L(2)-L(1);
%     %     axs(1).YLim = [(L(1)-Lgap*.05) (L(2)+Lgap*.05)];
%     axs(1).YLim = [(L(1)) (L(2))];
%     TickNums = linspace(L(1),L(2),minticks);
%     set(axs(1),'YTick',TickNums)
%     for l = 1:size(TickNums,2) 
%         axs(1).YTickLabel{l} = num2str(round(TickNums(l)*10^-(FirstNonEqualSigFig-2))*10^(FirstNonEqualSigFig-2)); 
%     end
% end

% % Update this later...
% if length(axs)>1
%     
%     L = [min(ydata2) maxs(ydata2)]; %
%     Lgap = L(2)-L(1);
%     if Lgap > 0
%         axs(2).YLim = [(L(1)-Lgap*.05) (L(2)+Lgap*.05)];
%         TickNums = linspace(L(1)-Lgap*.05,L(2)+Lgap*.05,minticks);
%         set(axs(2),'YTick',TickNums)
%         FirstNonEqualSigFig = log10(10.^floor(log10(abs(Lgap))));
%         if FirstNonEqualSigFig<-2 FirstNonEqualSigFig = -2; end
%         for l = 1:size(TickNums,2) axs(2).YTickLabel{l} = num2str(floor(TickNums(l)*10^-(FirstNonEqualSigFig-1))*10^(FirstNonEqualSigFig-1)); end
%     end
% end