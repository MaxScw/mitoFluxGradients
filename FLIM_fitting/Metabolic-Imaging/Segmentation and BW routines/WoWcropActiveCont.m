function [Wmask,acont,CoM,rad] = WoWcropActiveCont(im,NumWells)
% Takes in an image of a WoW (NADH OR FAD) dish and finds the mask of the 
% outer well to crop out, using active contours. The routine starts with a 
% disk of known radius, centered in the FOV center.
% clear all;
% im = tiffread2('C:\Users\Tim\Documents\Academic - Research\Data\2017-04-20 IllProfWells_test\2017-07-31 Shitty wells\s1_a1\sorted_sdts\Pos0_zBin_FADonly_IllProfed.tif');
% im = im.data;

% Versions:
% 2018-06-08: I updated this with the parameters most recently used in
%   MultiD_tiff_convert.m

if ~exist('NumWells')|NumWells==-1 NumWells = 9; end % Unless specified, assume 9-well dishes. Enter '16' if you're using 16-well dishes
if NumWells==9
    WellRad = 210;
elseif NumWells==16
    WellRad = 150;
else
    WellRad = NumWells; % If not 9 or 16, assume custom radius entered
end
xdim = size(im,2); ydim = size(im,1);
Gwell = fspecial('gaussian',[10 10],10);
gim = uint8(imfilter(InvertIm(im),Gwell,'same'));
Cx = round(xdim/2); Cy = round(ydim/2);
circ = Make_Circ_Mask(xdim,ydim,-1,-1,WellRad);
acont = activecontour(gim,circ,75);  %,'edge','ContractionBias',.2);
rad = sqrt(length(find(acont))./pi);
[x,y] = imCoM((double(acont)));
CoM = [x,y];

% NOTE: omit about 25 pixels around the edge because intensity
% attenuation happens 35 pixels in from the edge. Ie, avoid ~2/3 of
% the innacurate intensity area, but keep the 1/3 where attenuation
% isn't that bad so that we don't loose too much signal.
Wmask = Make_Circ_Mask(xdim,ydim,x,y,rad-25);


%% Optional plotting
% close all;
% imshow(im,[]); set(gcf,'position',[200 200 512 512])
% set(gcf,'paperPositionMode','auto')
% iptsetpref('ImshowBorder','tight');
%
% figure;
% imshow(gim,[]); set(gcf,'position',[200 200 512 512]);
% set(gcf,'paperPositionMode','auto')
% iptsetpref('ImshowBorder','tight');
%
% figure;
% imshow(circ,[]); set(gcf,'position',[200 200 512 512]);
% set(gcf,'paperPositionMode','auto')
% iptsetpref('ImshowBorder','tight');
%
% figure;
% imshow(a,[]); set(gcf,'position',[200 200 512 512]);
% set(gcf,'paperPositionMode','auto')
% iptsetpref('ImshowBorder','tight');
%
% [my,mx] = find(bwperim(Wmask));
% hold on;
% plot(mx,my,'r.')
