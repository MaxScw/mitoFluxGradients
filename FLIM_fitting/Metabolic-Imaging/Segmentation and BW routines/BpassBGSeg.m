function [masks singlemask num] = BpassBGSeg(im,lowfrq,hifrq,thresh,DilL,area_cut)
% Segmentation routine for performing a 1st-order segmentation using a
% bandpass filter and Otsu threshold.
% INPUTS:
% -im: input image
% -lowfrq: lower frequency cutoff (inverse of upper length scale)
% -hifrq: upper frequency cutoff (inverse of lower length scale)
% -thresh: lower pixel value cuttoff to threshold bpim (typically 1)
% -DilL: length of binary dilation after threshold

1;

if ~exist('lowfrq')|lowfrq==-1 lowfrq=2; end
if ~exist('hifrq')|hifrq==-1 hifrq=25; end
if ~exist('thresh')|thresh==-1 thresh=1; end
if ~exist('DilL')|DilL==-1 DilL=10; end
if ~exist('area_cut')|area_cut==-1 area_cut=[50 10^8]; end

bwSnipBool = 0; % Enable if you want to use bwSnipHighCurv.m 

bpim = gaussianbpf(im,lowfrq,hifrq);
bpmsk = zeros(size(im));
bpmsk(bpim>thresh)=1;
bpmsk = bwareafilt(boolean(bpmsk),area_cut);
bpmsk = imdilate(bpmsk,strel('disk',DilL));

% Sort into a cell
masks = {};
singlemask = bpmsk;
[L num] = bwlabel(bpmsk);
for b = 1:max(L(:))
    ind = find(L==b);
    [y,x] = find(L==b);
    tmp = zeros(size(im));
    tmp(ind) = 1;
    % Snip high curvature pinch points where embryos may be touching
    if bwSnipBool
        try
            snipmask = bwSnipHighCurv(tmp);
        catch
            snipmask = tmp;
        end
    else
        snipmask = tmp;
    end
    snipmask = bwareafilt(boolean(snipmask),area_cut);
%     close all
%     imshow(tmp); figure; imshow(snipmask)
    [L2 num2] = bwlabel(snipmask);
    for b2 = 1:max(L2(:))
        ind2 = find(L2==b2);
        tmp2 = zeros(size(im));
        tmp2(ind2) = 1;
        % Populate final masks cell
        masks = [masks tmp2];
    end
end

% Sort by x CoMs
num = size(masks,2);
XCoMs=[];
for b = 1:num
    [y,x] = find(masks{b});
    XCoMs(b) = mean(x);
end
[B,IX] = sort(XCoMs);
masks = masks(IX);

% Change to 1 to plot all these results
if 0
    close all; figure('position',[100 100 1000 800]);
    [y,x] = find(bpmskperim(bpmsk));
    subplot(2,2,1); imshow(im,[0 10]);
    subplot(2,2,2); imshow(bpim,[0 7]);
    subplot(2,2,3); imshow(bpmsk);
    subplot(2,2,4); imshow(im,[0 10]); hold on; plot(x,y,'b.')
end