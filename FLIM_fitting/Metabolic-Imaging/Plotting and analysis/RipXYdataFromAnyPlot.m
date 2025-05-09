function [x1, y1, z1] = RipXYdataFromAnyPlot(h)
% script gets the x and y data from all the current line
% plots in a figure (using gcf)

% if ~exist('h') h = gca; end
% dataObjs = get(h, 'Children'); %handles to low-level graphics objects in axes

if ~exist('h') h = gcf; end
dataObjs = findobj(h,'type','line'); %handles to low-level graphics objects in axes

% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object

x1 = get(dataObjs, 'XData');  %data from low-level grahics objects
y1 = get(dataObjs, 'YData');
z1 = get(dataObjs, 'ZData');