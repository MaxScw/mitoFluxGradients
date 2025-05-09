function HWPPowerCalMet(xlsfile)
% Input an excel file with only 4 measured calibration values.
% Protocol:
% 1) set Halfwave plate (HWP) to the angle that maximized the power through
% the objective. Calibrate that power to whatever your roof is (mine is
% 60mW at 750nm). Enter that as an angle with -1 in the same spreadsheet
% 2) Switch power meter to in front of the scanner and measure that power
% to get the transmission percentage between before-scanner and obj.
% 3) Set angle to 10, 18, 27, and 35 degrees, measure the power at the
% power meter (before scanner).
% Point this function to the excel spreadsheet, and it will interpolate the
% rest of the values based on the function I = Imax*sin(theta_max)^2. Will
% resave values to a new xls file.
% XLS columns = [angle,before scanner power]
% Put 'NADH' or 'FAD' in the xls filename.
% NOTE: 'Met' in title is because this is the version adapted to the new
% metabolic imaging scope

% if ~isempty(strfind(xlsfile,'NADH'))
%     trans = 0.095038009;
% elseif ~isempty(strfind(xlsfile,'FAD'))
%     trans = 0.072361473;
% else
%     error('Enter a proper channel, dumbass.')
% end

tab = xlsread(xlsfile);
LasOP{1} = -1;
LasOP{2} = tab(end,2);
tab(7:end,:)=[];
angs = tab(tab(:,1)>0,1)/1000; pows = tab(tab(:,1)>0,2);

Isin = @(par0,x) par0(1).*sin(2.*(x-par0(2))*pi/180).^2+par0(3);
% par0(1) = Imax;
% par0(2) = theta offset;
% Define initial parameters.  Assume amplitutde to be max-min of function,
% mux to be in the middle of the array, sigma to be 1/4 the array length,
% and the offset to be the mean of the function.
par0(1) = max(pows);
par0(2) = 0; % this can be optimized with center of mass.
par0(3) = 8; % Max extinction gives about 8mW at the objective

options = optimset('Display','off');   %changed from iter
par = lsqcurvefit(Isin,par0,angs,pows,[0 -45],[10000 45],options);
x = 0:50;
xfit = par(1).*sin(2.*(x-par(2))*pi/180).^2+par(3);
close all;
plot(angs,pows,'ro');
xlabel('angle'); ylabel('Power Before Objective')
hold on;
plot(x,xfit,'k');

outpows = [round(par(3)):300]';
% outangs = interp1(xfit,x,outpows);
outangs = (asin(sqrt((outpows-par(3))./par(1)))/pi*180/2+par(2));
plot(outangs,outpows,'g.')
spr = num2cell([outpows outangs.*1000]);
spr = [spr; LasOP];
dashes = strfind(xlsfile,'_');
xlswrite([xlsfile(1:dashes(end)-1) '_out.xlsx'],spr)
tb = cell2table(spr);
writetable(tb,[xlsfile(1:dashes(end)-1) '_out.txt'],'Delimiter','\t','WriteVariableNames',0);

