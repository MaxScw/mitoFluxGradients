function HWPPowerCal(xlsfile)
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
% rest of the values based on the function I = Imax*cos(theta_max)^2. Will
% resave values to a new xls file.
% XLS columns = [angle,before scanner power]
% Put 'NADH' or 'FAD' in the xls filename.

% if ~isempty(strfind(xlsfile,'NADH'))
%     trans = 0.095038009;
% elseif ~isempty(strfind(xlsfile,'FAD'))
%     trans = 0.072361473;
% else
%     error('Enter a proper channel, dumbass.')
% end

tab = xlsread(xlsfile);
angs = tab(tab(:,1)>0,1)/1000; pows = tab(tab(:,1)>0,2);
trans = tab(tab(:,1)<0,2)/max(pows);

Icos = @(par0,x) par0(1).*cos(2.*(x-par0(2))*pi/180).^2;
% par0(1) = Imax;
% par0(2) = theta offset;
% Define initial parameters.  Assume amplitutde to be max-min of function,
% mux to be in the middle of the array, sigma to be 1/4 the array length,
% and the offset to be the mean of the function.
par0(1) = max(pows);
par0(2) = 0; %this can be optimized with center of mass.

options = optimset('Display','off');   %changed from iter
par = lsqcurvefit(Icos,par0,angs,pows,[0 -45],[10000 45],options);
x = 0:55;
xfit = par(1).*cos(2.*(x-par(2))*pi/180).^2;
close all;
plot(angs,pows,'ro');
xlabel('angle'); ylabel('Objective Power')
hold on;
plot(x,xfit,'k');

outpows = [0 0.5 1 1.5 2:150]';
outpowsbeforesc = outpows./trans;
% outangs = interp1(xfit,x,outpowsbeforesc);
outangs = (acos(sqrt(outpowsbeforesc./par(1)))/pi*180/2+par(2));
plot(outangs,outpowsbeforesc,'g.')
spr = num2cell([outpows outpowsbeforesc outangs.*1000]);
dashes = strfind(xlsfile,'_');
xlswrite([xlsfile(1:dashes(end)-1) '_out.xlsx'],spr)

