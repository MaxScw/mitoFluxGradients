function snippedbw = bwSnipHighCurv(mask,MinPeakHeight,JointDiffsCut)
% Parameterize a perimeter with theta(s) contour coordinates. Look for high
% curvature regions and snip them off. Used for separating polar body form
% oocyte.

% Inputs:
% mask - binary mask to be analyzed for pinch points that should be snipped
% MinPeakHeight - peak height in units of rads/pixel step. 0.065 is a good value
% JointDiffsCut - High curvature points are identified below, then relative
%  cartesian and angular distance is calculated between pairs. We're
%  looking for points that are close, and close to anitparallel.
%  JointDiffsCut is a cutoff to exclude wiggles that aren't actually
%  antiparallel touching points between embryos. .6 is a good value for
%  512x512 images with 2X BH zoom.

% load('C:\Users\Tim\Documents\MATLAB\temp')
% mask = mask2;

% Versions:
% 2016-05-15: Refine snipping algorithm. Smooth binary masks. Also look for
% pinches that are surrounded by smooth regions (ie, peaks in theta that
% don't have other peaks close to them). Then further look for pairs of
% these peaks that are physically close to each other in xy.

warning('off','signal:findpeaks:largeMinPeakHeight')
if ~exist('MinPeakHeight')|MinPeakHeight==-1 MinPeakHeight=0.075; end
if ~exist('JointDiffsCut')|JointDiffsCut==-1 JointDiffsCut=.45; end
% MinPeakHeight=0.075; JointDiffsCut=.35; 


snippedbw = mask;
mask2 = imfill(mask,'holes');
G = fspecial('gaussian',[7 7],7);
mask2 = imfilter(mask2,G,'same');
% level = graythresh(double(mask2));
% % mask2 = im2bw(mask2,level);
% Make trace
skel = bwperim(mask2);

% Eliminate spurs and branch points
skel = bwmorph(skel,'spur');
skel = bwmorph(skel,'thin');
% Get xy coordinates from boundary trace. Coords in clockwise direction.
% Lowest y pixel is first pixel
B = bwboundaries(skel,'noholes');
% Pick largest region if a small region was created
blobsize = 0; biggest = 1;
for i = 1:length(B)
    if blobsize<size(B{i},1) 
        blobsize = size(B{i},1); biggest = i;
    end
end
coords = B{biggest};

% First find tangent vectors, then convert into theta.
tans = [diff(smooth(coords(:,1),20)) diff(smooth(coords(:,2),20))];
[theta, rho] = cart2pol(tans(:,1),tans(:,2));
[ThetAntiPar, rho] = cart2pol(-tans(:,1),-tans(:,2));
theta_unshifted = theta;
%  = theta+pi;
for i = 1:length(theta)-1
    if (theta(i+1)-theta(i))<-pi
        theta(i+1:end) = theta(i+1:end) + 2*pi;
    end
    if (theta(i+1)-theta(i))>pi
        theta(i+1:end) = theta(i+1:end) - 2*pi;
    end
    %     if (ThetAntiPar(i+1)-ThetAntiPar(i))<-pi
    %         ThetAntiPar(i+1:end) = ThetAntiPar(i+1:end) + 2*pi;
    %     end
    %     if (ThetAntiPar(i+1)-ThetAntiPar(i))>pi
    %         ThetAntiPar(i+1:end) = ThetAntiPar(i+1:end) - 2*pi;
    %     end
end
% plot(diff(theta))
% pks = findpeaks(-theta)

df = smooth(diff(theta),.05,'loess');
% Edges cause artifacts
df(1:10)=0; df(end-10:end)=0;

[pks,locs] = findpeaks(df,'MinPeakHeight',MinPeakHeight);
% Angles vary quickly, so do a 1st order interpolation to find the
% sub-pixel peaks.
for i=1:size(locs,1)
    [xm ym A]=crit_interp_p(df(locs(i)+(-2:2:2)),locs(i)+(-2:2:2));
    locsint(i) = xm;
    % get thetas for these peaks points
%     pkcoords(i,:) = [linterp(1:size(coords,1),coords(:,1),locsint(i)) linterp(1:size(coords,1),coords(:,2),locsint(i))];
    thetaint(i) = linterp(1:length(theta_unshifted),theta_unshifted,locsint(i));
    ThetAntiParint(i) = linterp(1:length(ThetAntiPar),ThetAntiPar,locsint(i));
end

% Get xy coords for peaks
pkcoords = [coords(locs,1) coords(locs,2)];

if size(locs,1)>1 % At least one pair of pinch points
    for i = 1:size(locs,1)
        diffs = sqrt((pkcoords(:,1)-pkcoords(i,1)).^2+(pkcoords(:,2)-pkcoords(i,2)).^2);
        ThetDiffs = abs(ThetAntiParint'-thetaint(i));
        diffs(i)=10^6; ThetDiffs(i)=10^6;
        JointDiffs = diffs./size(mask,1).*pi + ThetDiffs;
%         nearneighSpace(i,:) = [i find(diffs==min(diffs))];
%         nearneighThet(i,:) = [i find(ThetDiffs==min(ThetDiffs))];
        bestpair = find(JointDiffs==min(JointDiffs));
        bestpair = bestpair(1);
        nearneigh(i,:) = [i bestpair];
        JointDiffsFin(i) = JointDiffs(nearneigh(i,2));
    end
    
    % Filter pairs that aren't close enough in space and angle:
    nearneigh = nearneigh(JointDiffsFin<JointDiffsCut,:);
    if any(nearneigh)
        for i = 1:size(nearneigh,1)
            %                 imshow(snippedbw);
            x1 = pkcoords(nearneigh(i,1),2); y1 = pkcoords(nearneigh(i,1),1);
            x2 = pkcoords(nearneigh(i,2),2); y2 = pkcoords(nearneigh(i,2),1);
            xn = abs(x2-x1);
            yn = abs(y2-y1);
            
            % interpolate against axis with greater distance between points;
            % this guarantees statement is the under the first point!
            if (xn > yn)
                xc = x1 : sign(x2-x1) : x2;
                yc = round( interp1([x1 x2], [y1 y2], xc, 'linear') );
            else
                yc = y1 : sign(y2-y1) : y2;
                xc = round( interp1([y1 y2], [x1 x2], yc, 'linear') );
            end
            ind = sub2ind( size(mask), yc, xc );
            snippedbw(ind) = 0;
            % Make the cut 2-pixels thick to separate the new clusters
            ind = sub2ind( size(mask), yc, xc+1 );
            snippedbw(ind) = 0;
            
        end
    end
end
    
    