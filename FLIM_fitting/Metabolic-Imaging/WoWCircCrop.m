function [CropIm,crpind] = WoWCircCrop(im,Cx,Cy,r)
% Small function to identify a large circular region in the center of an
% image. Well-of-well (WoW) dishes have polystyrene autofluorescence around
% the well. This can be used to crop out this outer region with a circular
% mask of fixed radius.
% im - image
% Cx, Cy - center of disk
% r - disk radius
if ~exist('r')|r==-1 r = 215; end

xdim = size(im,2); ydim = size(im,1);
[Xm,Ym] = meshgrid((1:xdim)-Cx,(1:ydim)-Cy);
crpind = find(sqrt(Xm.^2+Ym.^2)>r);
CropIm = im;
CropIm(crpind) = 0;
