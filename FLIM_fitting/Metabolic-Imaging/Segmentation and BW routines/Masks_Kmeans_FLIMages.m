function [masks singlemask num] = Masks_Kmeans_FLIMages(im,FLim,Gblur,nClusts,KeepClusts,area_cut,DilPix,circ_cut,PrevCoMs,crpind)
% Take in an image, threshold it, then blur, then thress. Then use
% consecutive erodes to separate blobs, then use 'watershed' process to
% find final masks
% PrevCoMs - if previous egg CoMs are included, this runs a radius filter,
%   blocking out all pixels not within a radius of any of the eggs. Radius
%   is sqrt(area/pi).

% clear all;
% area_cut = [1000 100000];
% Gblur = 6;
% im = imread('Z:\Lab\Tim\2016-04-23 Old Young Embs\2016-04-29 Repeat\s1_Day1\sorted_sdts\IntTiffs_Pos3_s1_Day1\fr00120.tif');
% FLim = imread('Z:\Lab\Tim\2016-04-23 Old Young Embs\2016-04-29 Repeat\s1_Day1\sorted_sdts\FLIMageTiffs_Pos3_s1_Day1\fr00120.tif');
% nClusts = 4;
% KeepClusts = [2 3 4];

% Versions:
% 2018-10-06: Repurposed as initial pass before Ilastik segmentation. Also
% altered to exclude cropped pixels from kmeans grouping.
% 2015-07-16: Was trying to segment the nucleus with a secondary kmeans
% after embryo-vs-bg kmeans. Didn't work well. Instead just do a 3-color
% kmeans in the beginning to find the nuclei.
% 2015-07-14: Adapted from 'EllMask_blur_fillholes', but use gauss blur and
% k-means clustering to get masks instead.

bwSnipBool = 0; % Enable if you want to use bwSnipHighCurv.m 

% For NADH*FAD images, 3 tends to get: nucleus and background, dim cytosol,
% bright cytosol
xdim = size(im,2);
ydim = size(im,1);
if ~exist('nClusts')|nClusts==-1 nClusts = 3; end % Number of layers for k-means. 
if ~exist('KeepClusts')|KeepClusts==-1 KeepClusts = [2 3]; end % Number of layers for k-means. 
if ~exist('circ_cut')|circ_cut==-1 circ_cut = 10; end % Unless specified, assume amorphous shapes are OK
if ~exist('bpcut') bpcut = [3 70]; end
if ~exist('FLim')|FLim==-1
    FLim = im;
    FLimBool = 0;
end
if ~exist('crpind')|crpind==-1 crpind = []; end
G = fspecial('gaussian',[Gblur Gblur],Gblur);


% If NADH and FAD present, combine into one image by multiplying together
if xdim==2*ydim
    FADBool = 1;
    NADHBool = 1;
    NADim = double(im(:,1:xdim/2)); NADFLim = double(FLim(:,1:xdim/2));
    FADim = double(im(:,xdim/2+1:end)); FADFLim = double(FLim(:,xdim/2+1:end));
    xdim = size(NADim,2); ydim = size(NADim,1);
    im = (NADim./mean(NADim(:))).*(FADim./mean(FADim(:)));
    FLim = (NADFLim./mean(NADFLim(:))).*(FADFLim./mean(FADFLim(:)));

%     im = (NADim./mean(NADim(:)))+(FADim./mean(FADim(:)));
%     FLim = (NADFLim./mean(NADFLim(:)))+(FADFLim./mean(FADFLim(:)));
end
% Else run kmeans on whatever channel is present

% Optional radial filter if PrevCoMs is passed to function
if exist('PrevCoMs')
    RadCutMask = zeros(ydim,xdim);
    rad = sqrt(area_cut(2)/pi);
    for i = 1:size(PrevCoMs,1)
        [X,Y] = meshgrid((-PrevCoMs(i,1)+1):(-PrevCoMs(i,1)+xdim),(-PrevCoMs(i,2)+1):(-PrevCoMs(i,2)+ydim));
        Z = sqrt(X.^2 + Y.^2);
        Z(Z>rad) = 0; Z(find(Z))=1;
        RadCutMask = RadCutMask|Z;
    end
    ind = find(~RadCutMask);
    im(ind) = 0;
end

% k-means
gim = imfilter(im,G,'same');
kmint = double(gim);
kmint(crpind) = 0;
kmint = double(reshape(kmint,xdim*ydim,1));
kmint = (kmint-min(kmint))./max(kmint);

% Haven't used FLIMages for a long time. Comment out for now and just use
% ints.
% if FLimBool    
%     gFLim = imfilter(FLim,G,'same');
%     kmFLIM = double(gFLim);
%     kmFLIM = double(reshape(kmFLIM,xdim*ydim,1));
%     kmFLIM = (kmFLIM-min(kmFLIM))./max(kmFLIM);
%     kmarr = [kmint kmFLIM];
% else
%     kmarr = [kmint];
% end
% % Ocropped regions that are 0.
% kmarr(kmarr==0)=[];

[cluster_idx, cluster_center] = kmeans(kmint,nClusts,'distance','sqEuclidean', ...
    'Replicates',3);
[tmp,srtclusts] = sort(cluster_center);
keepinds = [];
cluster_idx_imdims = cluster_idx;
% cluster_idx_imdims = zeros(size(kmint));
% cluster_idx_imdims(NonCrp) = cluster_idx;
keeparr = zeros(size(cluster_idx_imdims,1),1);
for i = KeepClusts
    keeparr(cluster_idx_imdims==srtclusts(i))=1;
end
bw = reshape(keeparr,ydim,xdim); % make mask of pixels from cluster 2
%     imshow(bw,[])
%     [y,x] = find(bwperim(bw)); imshow(im,[0 500]); hold on; plot(x,y,'.r')
bw = bwareafilt(boolean(bw),area_cut);

% Sort into a cell
masks = {};
singlemask = bw;
[L num] = bwlabel(bw);
for b = 1:max(L(:))
    ind = find(L==b);
    [y,x] = find(L==b);
    tmp = zeros(ydim,xdim);
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
        tmp2 = zeros(ydim,xdim);
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
