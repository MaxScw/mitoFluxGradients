function BG_ints_plot(acqpath,pos)
% Short script to plot background intensities to check if they drift,
% either due to laser power drift or z drift.
 
% % TEST in script mode:
% clear all;
% acqpath = 'Z:\Lab\Marta\Flow Chamber Tests\2019-04-10_Flowtest_FCCP2\s1_a1';
% pos = 1;

if acqpath(end)~='\' acqpath = [acqpath '\']; end
if ~exist('pos')|pos==-1 pos = 1; end
ChLabs = {'NADH','FAD','UserChan'};
startupTim

close all;
load([acqpath 'BG_vals.mat'])
% BGvals matrix dims: [ValType, framenum, posnum, ch];
%    -> ValTypes: 1=intensity, 2=std, 3=NumOfPix, 4=timestamp

% Select appropriate time units. If more that 90 min, use hours
BGtimestps = squeeze(BGvals(4,:,:,:));
ExpT = (max(BGtimestps(BGtimestps>0))-min(BGtimestps(BGtimestps>0)))*86400/60;

% 60 gives minutes, 3600 for hours, 86400 for days
if ExpT < 90
    tlab = 'Time (min)';
    timeunit = 60;
else
    tlab = 'Time (h)';
    timeunit = 3600;
end

% % Optional (uncomment), look at masks
% close all; figure('color',[.5 .5 .5])
% for i = 1:size(BGmasks,1)
%     imshow(BGmasks{i},[]); pause(.2);
% end

BadFrs = find(BGvals(4,:,pos,1)==0);
T = squeeze(BGvals(4,:,pos,1));
T(BadFrs)=[]; T = (T-min(T))*86400/timeunit;
Ys = squeeze(BGvals(1,:,pos,:)); Ys(BadFrs) = [];

Ny = squeeze(BGvals(1,:,pos,1)); Ny(BadFrs) = [];
Fy = squeeze(BGvals(1,:,pos,2)); Fy(BadFrs) = [];

Nstd = squeeze(BGvals(2,:,pos,1))./squeeze(sqrt(BGvals(3,:,pos,1))); Nstd(BadFrs) = [];
Fstd = squeeze(BGvals(2,:,pos,2))./squeeze(sqrt(BGvals(3,:,pos,2))); Fstd(BadFrs) = [];

if ~isempty(Ny)&~isempty(Fy)
    h=figure('position',[462   500   658   448])
    [ax,h1,h2] =plotyy(T,Ny,T,Fy);
    xlabel(tlab);
    hold(ax(1),'on'); hold(ax(2),'on');
    % NADH
    e1 = errorbar(ax(1),T,Ny,Nstd);
    % FAD
    e2 = errorbar(ax(2),T,Fy,Fstd);
    ax(1).YLabel.String = 'NADH BG int'; ax(2).YLabel.String = 'FAD BG int';
    
    ax(1).YLim = [min(min(Ny)) max(max(Ny))];
    ax(2).YLim = [min(min(Fy)) max(max(Fy))];
    ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
    e1.Color = [0 0 1]; e2.Color = [0 1 0];
    set(gca,'position',[0.1459      0.12946      0.71277      0.79554])
    GoodPlotyyTicks(ax);
    saveas(gcf,[acqpath 'BG_masks_ints_Pos' num2str(pos-1) '.fig'])
    close(h)
else
    Uy = squeeze(BGvals(1,:,pos,3)); Uy(BadFrs) = [];
    Ustd = squeeze(BGvals(2,:,pos,3))./squeeze(sqrt(BGvals(3,:,pos,3))); Ustd(BadFrs) = [];
    % Not worth coding/fixing right now. Fix if needed.
%     for ch = 1:3
%         if ~isempty(squeeze(BGvals(1,:,pos,ch)))
%             h=figure('position',[462   500   658   448])
%             [ax,h1,h2] =plotyy(T,Ny,T,Fy);
%             xlabel(tlab);
%             hold(ax(1),'on'); hold(ax(2),'on');
%             % NADH
%             e1 = errorbar(ax(1),T,Ny,Nstd);
%             % FAD
%             e2 = errorbar(ax(2),T,Fy,Fstd);
%             ax(1).YLabel.String = 'NADH BG int'; ax(2).YLabel.String = 'FAD BG int';
%             
%             ax(1).YLim = [min(min(Ny)) max(max(Ny))];
%             ax(2).YLim = [min(min(Fy)) max(max(Fy))];
%             ax(1).YColor = [0 0 1]; ax(2).YColor = [0 1 0];
%             e1.Color = [0 0 1]; e2.Color = [0 1 0];
%             set(gca,'position',[0.1459      0.12946      0.71277      0.79554])
%             GoodPlotyyTicks(ax);
%             saveas(gcf,[acqpath 'BG_masks_ints_Pos' num2str(pos-1) '_' ChLabs{ch} '.fig'])
%             close(h)
%         end
%     end
end




