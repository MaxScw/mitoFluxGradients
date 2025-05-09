function coords = Skel2NeighborCoords(mask)
% Takes in a skeletenized (1-pixel-thin) binary image, often produced using
% 'bwperim'. This function sorts the pixels into coordinates that are
% ordered by sequential neighbors, walking along the contour.

% clear all; load('C:\Users\Tim\Desktop\test\test.mat')

% Versions:
% 2016-05-09: Use Matlab's much quicker 'bwboundaries' to get orderedcoords

skel = bwperim(imfill(mask,'holes'));

% Eliminate spurs and branch points
skel = bwmorph(skel,'spur');
skel = bwmorph(skel,'thin');

B = bwboundaries(skel,'noholes');
coords = B{1};

% [b,a] = find(bwmorph(skel,'branchpoints'));
% skel(ind) = 0;



