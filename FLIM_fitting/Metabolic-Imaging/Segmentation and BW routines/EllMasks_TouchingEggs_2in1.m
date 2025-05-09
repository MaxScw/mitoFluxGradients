function [acs matches] = EllMasks_TouchingEggs_test(im,egg,ma)
% Takes an image with a number of eggs touching. Returns a structure of
% masks of each egg.
% -Performs successive erodes to separate eggs from each other and find
%  centers of masses
% -Rotates an egg template image and convolutes locally around each center
%  of mass to find best orientation and overlap
% -Uses an active contour to get a tighter boundary around each specific
% egg
% im - the image
% egg - template of the egg
% ma - mask of the egg template

% clear all
% st = tiffread2('C:\Users\Tim\Documents\Academic - Research\Data\test\Pos0_IntenStk_s1NADH.tif'); im = st(1).data;
% ma = imread('C:\Users\Tim\Documents\Academic - Research\Data\test\mask.tif'); map = imread('C:\Users\Tim\Documents\Academic - Research\Data\test\mask_per.tif');
% egg = imread('C:\Users\Tim\Documents\Academic - Research\Data\test\eggave.tif');

map = bwperim(ma);

%% Erode main image
xdim = size(im,2); ydim = size(im,1);
no = wiener2(im,[5 5]);
level = graythresh(im);
bw = im2bw(im,level);
f = imfill(bw,'holes');

% Find connected clusters of pixels and perimeters
er = imerode(f,strel('disk',10));

% % Iterative erode. Unnecessary, I think. Just once more and it seems to do
% % very well and counting the number of eggs for you.
% [L, num] = bwlabel(er);

% i=1;
% while num<numeggs
%     i = i+1;
%     er = imerode(er,strel('disk',10));
%     [L, num] = bwlabel(er);
%     subplot(2,2,3); imshow(er);
%     pause(.1)
% end

er = imerode(er,strel('disk',10));
[L, num] = bwlabel(er);

% % Display
% close all;
% figure('units','normalized','position',[.1 .1 .8 .8]);
% subplot(2,2,1); imshow(im,[0 10])
% subplot(2,2,2); imshowt(f)

%% Rotate mask by 10deg intervals and convolute with the whole image, find
% local peaks to compare with CoMs
peaks=zeros(num,3,36);
for j = 1:36
    close all;
    rots(j) = (j-1)*10;
    eggr = double(imrotate(egg,rots(j))); % 10 deg rotation increments
    %         G = fspecial('gaussian',[10 10],8);
    %         eggr = imfilter(eggr,G,'same');
    eggr = eggr - mean(mean(eggr));
    
    C = conv2(im,eggr,'same');
    % Find peaks iteratively
    rad = min(size(egg,1),size(egg,2))*.5;
    tmppeaks = findpeaksT(C,num,rad);
    peaks(:,:,j) = tmppeaks;
    j
end

%% Find best convolution match and rotation for each egg
% Compare peaks to CoM's of each blob/egg. Points must be within 1/10 the
% radius to qualify. Then find the point with the highest score
best_matches=zeros(num,4);
% subplot(2,2,4)
ihm = imshow(im,[0 10],'border','tight'); hold on;


for i = 1:num
    clear scores
    [y,x] = find(L==i);
    CoMs(i,1) = mean(x);
    CoMs(i,2) = mean(y);
    % Calculate distances in a 2d matrix of dims [peaknum,rotnum]
    dists = squeeze(sqrt((peaks(:,1,:)-CoMs(i,1,:)).^2+(peaks(:,2,:)-CoMs(i,2,:)).^2));
    [peaknums,rotnums] = find(dists<rad/10);
    for j = 1:size(peaknums,1)
        scores(j) = peaks(peaknums(j),3,rotnums(j));
    end
    bestind = find(scores==max(scores));
    best_matches(i,:) = [peaks(peaknums(bestind),:,rotnums(bestind)) (rotnums(bestind)-1)*10];
    % columns: [x y score rot]
    
    % Plot best matches over the image
    mar = double(imrotate(ma,(rotnums(bestind)-1)*10));
    %     mar = imerode(mar,strel('disk',3));
    mapr = imrotate(map,(rotnums(bestind)-1)*10);
    [my,mx] = find(mapr);
    mx = mx + best_matches(i,1) - round(size(mapr,2)/2);
    my = my + best_matches(i,2) - round(size(mapr,2)/2);
    plot(mx,my,'r.','markersize',3)
    matches{i} = [mx my];
    
    % Make masks for active contours
    ind1 = (best_matches(i,2)-round(size(mar,1)/2)+1):(best_matches(i,2)+round(size(mar,1)/2)-1);
    ind2 = (best_matches(i,1)-round(size(mar,2)/2)+1):(best_matches(i,1)+round(size(mar,2)/2)-1);
    acmask{i} = zeros(size(im,1),size(im,2));
    acmask{i}(ind1,ind2) = mar;
end

%% Lastly, make it a tight boundary with an active contour
acmask0 = acmask;
for i = 1:num
    overlap = zeros(size(im,1),size(im,2));
    for k = 1:num
        if i~=k
            acmask{i}(find(acmask0{i}==acmask0{k}))=0;
            %             acmask{i}(find(acmask0{i}==imerode(acmask0{k},strel('disk',2))))=0;
        end
    end
    ind = find(overlap);
    acmask{i}(ind)=0;
    pause(0.5)
    
    acs{i} = activecontour(no,imerode(acmask{i},strel('disk',2)),200,'Chan-Vese',2);
    acp = bwperim(acs{i});
    [y,x] = find(acp);
    plot(x,y,'g.','markersize',3)
    set(gcf,'PaperPositionMode','auto')
    
end

