% function errorbarTS(ax,x,y,yerr)
% I had to write a custom error bar routine because stupid matlab doesn't
% all for partial transparency of the 'errorbar' object type

% Base on https://stackoverflow.com/questions/43805596/draw-error-bars-in-matlab
% ... didn't complete, too much labor. Stupid Matlab.

cl = 1;
figure;
ax = gca;
x = T(:,cl);
y = dat(:,cl);
yerr = daterr(:,cl);
y_limits = [y-yerr./2 y+yerr./2]

for i = 1:length(x)
ln = line([x(i) x(i)], [y_limits(i,1) y_limits(i,2)]);
ln
end