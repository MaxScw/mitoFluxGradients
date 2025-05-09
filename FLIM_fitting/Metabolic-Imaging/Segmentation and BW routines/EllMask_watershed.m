
function [masks singlemask num] = EllMask_watershed(NADim,FADim,area_cut,level_fact)
% Take in an image, threshold it, then blur, then thress. Then use
% consecutive erodes to separate blobs, then use 'watershed' process to
% find final masks

% clear all;
% level_fact=1;
% area_cut = [100 50000];

% im = imread('C:\Users\Tim\Documents\Academic - Research\Data\2014-05-09 Cumulus cells first\testNADHim1.tif');
% load('C:\Users\Tim\Documents\Academic - Research\Data\s2m1_VC2290\tempdat');
% im = NADim;

% Versions:
% -Working for mouse oocytes (big blobs), but not for C elegans. Also, 
% sometimes FAD is bright, NADH is dim or vice-versa. Find all eggs by 
% doing thresh in both channels, then combine them.

xdim = size(NADim,2); ydim = size(NADim,1);
if ~exist('feat_diam') egg_diam = ydim/3; end
if ~exist('circ_cut') circ_cut = 10; end % Unless specified, assume amorphous shapes are OK

im = NADim;
% im = histeq(NADim);
% im = uint8((im-min(im(:)))./max(im(:)).*255);
% K = wiener2(im,[5 5]);
K = bpassTS(NADim,3,70);
level = graythresh(K);
bw1 = im2bw(K,level*level_fact);
G = fspecial('gaussian',[15 15],10);
gim = imfilter(double(bw1),G,'same');

level = graythresh(gim*level_fact);
bw = im2bw(gim,level*level_fact);

% bw = imfill(bw,'holes'); % close all;imshowt(bw);figure;imshowt(im)
Nbwe = bw; 

im = FADim;
% im = histeq(FADim); 
% im = uint8((im-min(im(:)))./max(im(:)).*255);
% im = histeq(double(FADim)); im = uint8((im-min(im(:)))./max(im(:)).*255);
if ~exist('feat_diam') egg_diam = ydim/3; end
if ~exist('circ_cut') circ_cut = 10; end % Unless specified, assume amorphous shapes are OK
% K = wiener2(im,[5 5]);
K = bpassTS(FADim,3,70);
level = graythresh(K);
bw1 = im2bw(K,level*level_fact);
G = fspecial('gaussian',[15 15],10);
gim = imfilter(double(bw1),G,'same');

level = graythresh(gim*level_fact);
bw = im2bw(gim,level*level_fact);

% bw = imfill(bw,'holes'); % close all;imshowt(bw);figure;imshowt(im)
Fbwe = bw;

% Merge masks for both NADH and FAD! Sometimes brightness disparity causes
% algorithm to miss eggs in one channel or other
bw_0 = Nbwe|Fbwe;
% bw_0 = bw;

% Do erodes to separate touching eggs. Erode each cluster separately
% until cluster dissapears
bwe = bw_0;
% Put each blob into a structure element
StableBlobs = [];
[L num] = bwlabel(bwe);
for b = 1:max(L(:))
    % Quick small area filter
    ind = find(L==b);
    [y,x] = find(L==b);
    if size(ind,1)<area_cut(1)
        bwe(ind)=0;
    else
        Blob.mask = zeros(ydim,xdim);
        Blob.finmask = zeros(ydim,xdim);
        Blob.mask(ind) = 1;
        Blob.CoM =[mean(x) mean(y)];
        Blob.Erodes = 0; 
        % Keeps track of total times blobs have been erorded. Updated when
        % blobs split. Used to dilate the blobs after separation
        StableBlobs =[StableBlobs Blob];
%         imshow(Blob.mask)
    end
    
end
NewBlobs = [];

% Do erodes independently for each blob. If a blob splits, start over for
% both of the new blobs
Remove = [];
BreakBool=0;
bnum = 1;
while bnum <= size(StableBlobs,2)
%     imshow(StableBlobs(bnum).mask)
    erblob = StableBlobs(bnum).mask;
    eri = 0;
    num = 1;
    while num>0
        erblob = imerode(erblob,strel('disk',2));
        eri = eri + 1;
%         imshow(erblob)
%         pause(.05)
        [L num] = bwlabel(erblob);
        
        if num>1
            % Quick small area filter. Only keep blobs that are larger that
            % area cutoff
            CutBlobs = 0;
            for i = 1:num
                inds{i} = find(L==i);
                areas(i) = size(inds{i},1);
                if areas(i)<area_cut(1)
                    erblob(inds{i})=0;
                    CutBlobs = CutBlobs + 1;
                end
            end
            % If all blobs in this erode fell below area cut, use the largest
            % one:
            if CutBlobs==num
                ind = inds{areas==max(areas)};
                erblob(ind(1))=1;
            end
            [L num] = bwlabel(erblob);
        end
        
        % If blob didn't split
        if num==1
            
        elseif num>1
            % If it did split into two blobs that are larger than area cut,
            % replace current blob with first child-blob
            % and send other child-blobs to the end of the StableBlobs
            % structure
            TotEris = StableBlobs(bnum).Erodes + eri;
            ind = find(L==1);
            [y,x] = find(L==1);
            StableBlobs(bnum).mask = zeros(ydim,xdim);
            StableBlobs(bnum).finmask = zeros(ydim,xdim);
            StableBlobs(bnum).Erodes = TotEris;
            StableBlobs(bnum).mask(ind) = 1;
            StableBlobs(bnum).CoM =[mean(x) mean(y)];
            for i=2:max(L(:))
                ind = find(L==i);
                [y,x] = find(L==i);
                Blob.mask = zeros(ydim,xdim);
                Blob.finmask = zeros(ydim,xdim);
                Blob.mask(ind) = 1;
                Blob.Erodes = TotEris;
                Blob.CoM =[mean(x) mean(y)];
                StableBlobs =[StableBlobs Blob];
            end
            bnum = bnum - 1;
            break;
        end
    end
    bnum = bnum + 1;    
end
% Now dilate the remaining blobs to compensate for the erosions at the time
% of split:
% imshow(NADim+FADim,[]); hold on;
for b = 1:size(StableBlobs,2)
    im = StableBlobs(b).mask;
    for i = 1:round(StableBlobs(b).Erodes*2/3);
        im=imdilate(im,strel('disk',1));
    end
    StableBlobs(b).mask = im;
    [y,x] = find(bwperim(im));
%     plot(x,y,'.r')
end

%% Watershed segment separating
D = bwdist(~bwe,'chessboard');
D = -D;
% subplot(2,2,2) 
% imshow(D,[])
D(~bwe) = -Inf;
L = watershed(D);
L = L-1;
% figure;imshow(L,[])

% Match all segments to the eroded mask that they have the most overlap.
for i = 1:max(L(:))
    for j = 1:size(StableBlobs,2)
        waterblob = zeros(ydim,xdim);
        waterblob(find(L==i)) = 1;
        overlap(j) = length(find(waterblob & StableBlobs(j).mask));
        [y,x] = find(waterblob);
        dists(j) = sqrt((mean(x)-StableBlobs(j).CoM(1))^2+(mean(y)-StableBlobs(j).CoM(2))^2);
    end
    % If eroded blob is too small, a watershed might not have any overlap.
    % In this case, find the closest CoM to match it up
    if isempty(find(overlap))
        blobmatch = find(dists==min(dists));
    else
        blobmatch = find(overlap==max(overlap));
    end
    StableBlobs(blobmatch).finmask = StableBlobs(blobmatch).finmask + waterblob;
end

singlemask = zeros(ydim,xdim);
num = size(StableBlobs,2);

for b = 1:num
    masks{b} = bwmorph(StableBlobs(b).finmask,'bridge');
    [y,x] = find(masks{b});
    XCoMs(b) = mean(x);
    singlemask = singlemask + masks{b};
end

[B,IX] = sort(XCoMs);
StableBlobs = StableBlobs(IX);
masks = masks(IX);

% %% Plot masks over image for diagnostic. (uncomment)
% figure; imshow(NADim+FADim,[]); hold on
% for b = 1:size(masks,2)
%     [y,x] = find(bwperim(masks{b}));
%     plot(x,y,'.','color',[b/size(StableBlobs,2) b/size(StableBlobs,2) 1-b/size(StableBlobs,2)]);
%     plot(StableBlobs(b).CoM(1),StableBlobs(b).CoM(2),'x','color',[b/size(StableBlobs,2) b/size(StableBlobs,2) 1-b/size(StableBlobs,2)]);
% end
% disp('')

% Previous version had circularity and whatnot, but deleted here
