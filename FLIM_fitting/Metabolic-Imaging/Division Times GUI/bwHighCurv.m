function highcurvpnts = bwHighCurv(mask,MinPeakHeight,JointDiffsCut,smo)
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
%  antiparallel touching points between embryos. 
%  0.6 is a good value for 512x512 images with 2X BH zoom.
% smo - parameter for smoothing differentiation, default to 15 seems good
% load('C:\Users\Tim\Documents\MATLAB\temp')
% mask = mask2;

% Versions:
% 2016-05-15: Refine snipping algorithm. Smooth binary masks. Also look for
% pinches that are surrounded by smooth regions (ie, peaks in theta that
% don't have other peaks close to them). Then further look for pairs of
% these peaks that are physically close to each other in xy.

warning('off','signal:findpeaks:largeMinPeakHeight')
if ~exist('MinPeakHeight')|MinPeakHeight==-1 MinPeakHeight=0.05; end
if ~exist('JointDiffsCut')|JointDiffsCut==-1 JointDiffsCut=.45; end
if ~exist('smo')|smo==-1 smo=11; end
%MinPeakHeight=0.075; JointDiffsCut=.35; 


snippedbw = mask;
mask2 = imfill(mask,'holes');
G = fspecial('gaussian',[7 7],7);
mask2 = imfilter(mask2,G,'same');
level = graythresh(mask2);
mask2 = im2bw(mask2,level);
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
tans = [diff(smooth(coords(:,1),smo)) diff(smooth(coords(:,2),smo))];
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
    thetaint(i) = interp1(1:length(theta_unshifted),theta_unshifted,locsint(i),'linear');
    ThetAntiParint(i) = interp1(1:length(ThetAntiPar),ThetAntiPar,locsint(i));
end



% Get xy coords for peaks
pkcoords = [coords(locs,1) coords(locs,2)];
%highcurvpnts=size(pks,1);
%highcurvpnts=locs;
if ~isempty(locs)
    highcurvpnts=[size(locs,1),mean(diff(locs))];
else
    highcurvpnts=[size(locs,1),0];
end
%plot(coords(:,1),coords(:,2))
%hold on;
%plot(coords(locs,1),coords(locs,2),'ro')
%hold off;

    