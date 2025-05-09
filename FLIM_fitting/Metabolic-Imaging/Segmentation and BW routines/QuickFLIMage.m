function flimage = QuickFLIMage(flim,Bin,thresh)
% Calculates a quick FLIM image by taking the mean lifetime for each pixel
% (BH refer to this as the 1st moment). Also allows a pre-Binning off
% photons from neighboring pixels to increase S:N for each pixel, same way
% BH does in SPCimage.
% INPUTS:
% -FLIM: the full imported FLIM file, typically [256,512,512] dimensions
% -Bin: number of neighboring pixels to add to current pixel
% -thresh: minimum number of photons (post-Bin) for calculating a mean.
%   Otherwise pixel has mean lifetime of 0.
% OUTPUTS:
% FLIMage: mean lifetime image. Pixel values are on same scale as input
%   nanotime, typically 256 Bins. Convert to real time units outside this
%   function using the TAC information, if desired.

% % TEST:
% clear all;
% load('tmp.mat')
% flim = double(nadch1);
% Bin = 1;
% thresh = 20;

% BIN photons (blur) to improve decays. 
flimbin = zeros(size(flim));
for px1 = 1:size(flim,2)
    for px2 = 1:size(flim,3)
        % adjusted in case pixel is at the edge of the image
        px1adj = [max([1 (px1-Bin)]) min([size(flim,2) (px1+Bin)])];
        px2adj = [max([1 (px2-Bin)]) min([size(flim,2) (px2+Bin)])];
        flimbin(:,px1,px2) = sum(sum(flim(:,px1adj(1):px1adj(2),px2adj(1):px2adj(2)),2),3);
    end
end

% Thresh
flimsum = squeeze(sum(flimbin,1));
[y,x] = find(flimsum>=thresh);

flimage = zeros(size(flim,2),size(flim,3));

% Calculate a FLIMage
for i = 1:length(x)
    ind = find(flimbin(:,y(i),x(i)));
    % weight bins by number of photons in them. Ie, find center of mass of dist
    flimage(y(i),x(i)) = sum(ind.*flimbin(ind,y(i),x(i)))/sum(flimbin(ind,y(i),x(i)));
%     plot(flimbin(:,y(i),x(i)))
%     [wtnanmean(ind,flimbin(ind,y(i),x(i))) wtnanmean(ind,flimbin(ind,y(i),x(i)))/256*10]
%     1;
end
% imshow(flimage,[])
