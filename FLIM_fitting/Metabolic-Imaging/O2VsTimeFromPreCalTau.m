function [O2T, O2] = O2VsTimeFromPreCalTau(O2WritePath,DropTime,RestTime,O2Cons,DecTau)
% Adapted from O2FromDecayTau to accept known O2 decay tau value, obtained
% from a previous O2 calibration run with same geometry as the experimental
% dish with embryos in it.

% Restore: O2 typically restored to original concentration at a later time.
% Include this info.

% Inputs:
% -O2WritePath: Path to write O2_vs_time data. Make this the same path as
%   the path of the corresponding O2drop acquisition you took.
% -DecTau: Tau measured by calibration measurements of Ruthenium dye.
%   Unites are in hours.
% -DropTime: O2 drop start finish times. [t0 t1] -> DateString formats
%   AT LEAST DROP RANGE IS REQUIRED FOR THIS VERSION
% -RestTime: O2 restore start finish times. [t0' t1']
% -O2Cons: Initial O2 concentration and drop concentration. [O2_i O2_f]

% Versions:
% 2016-12-28: Simplified Drop and Restore input times. Instead of ranges,
%  just enter the times when the Drop and Restore tank switches happened.

% NOTE: all times in absolute time units

% clear all
% DecTau = 4.9185/60;
% O2WritePath = 'Z:\Lab\Emily\2016-11-10 Emily First O2 Drop on Embryos\2016-12-02 Repeat\s1_a1_O2drop';
% DropTime = '2016-12-02 16:53:00';
% RestTime = '2016-12-02 18:05:00';
% O2Cons = [5 0];

% timeunit = 3600; % 60 gives minutes, 3600 for hours, 86400 for days
% switch timeunit
%     case 1;
%         tlab = 'Time (s)';
%     case 60;
%         tlab = 'Time (min)';
%     case 3600;
%         tlab = 'Time (h)';
%     case 86400;
%         tlab = 'Time (days)';
% end

% Versions:

if ~exist('DecTau')|DecTau==-1 DecTau = 4.9185/60; end
if ~exist('O2Cons')|O2Cons==-1 O2Cons = [5 0]; end


%% Convert times to datenum
Dsdt = dir([O2WritePath '\*.sdt']);
EndTime = Dsdt(end).datenum;
if exist('DropTime')
    DropTime = datenum(DropTime)+1/24/60;
    % Add a minute to everything because it takes about a minute to flush the
    % chamber with gas.
else
    error('Drop range must be specified for this routine, oaf!')
end
% And if not specified, assume no restore occurs. Set st/end both to end 
if exist('RestTime')
    RestTime = datenum(RestTime)+1/24/60;
else
    RestTime = EndTime;
end

%% Construct O2 vs time series with a simple exp equation


% Get O2 values for O2 drop range
DropT = 0:.005:(RestTime-DropTime)*24;
DropO2 = (O2Cons(1)-O2Cons(2))*exp(-DropT./DecTau)+O2Cons(2);
DropT = (DropT)/24+DropTime(1); % Back to datenum units

% Get O2 values for O2 restore range
RestT = 0:.005:(EndTime-RestTime)*24;
RestO2 = O2Cons(1)-(O2Cons(1)-O2Cons(2))*exp(-RestT./DecTau);
RestT = (RestT)/24+RestTime(1); % Back to datenum units

%% Stitch full O2 vs exp time (relative to 1st frame save time), and save
O2T = [DropT'; RestT'];
O2 = [DropO2'; RestO2'];

save([O2WritePath '\O2Info.mat'],'O2T','O2','DropTime','RestTime','EndTime','O2Cons','DecTau');

%% Plot O2 vs time and save
figure('position',[800 200 500 350]);
% figure('units','normalized','position',[.1 .1 .35 .3]);
plot((O2T-O2T(1))*24,O2);
set(gca,'fontsize',10)
xlabel('time (h)','fontsize',12);ylabel('O2 conc (%)','fontsize',12);
set(gcf,'paperPositionmode','auto')
saveas(gcf,[O2WritePath '\O2VsT_vals.jpg'])
saveas(gcf,[O2WritePath '\O2VsT_vals.fig'])

