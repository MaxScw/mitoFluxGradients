function [matchcoords acmasks numeggs] = EllMasks_TouchingEggs_ConvMatches(im)
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
load('C:\Users\Tim\Documents\Academic - Research\Data\tempdat'); 
% im = NADim;
map = bwperim(ma);

%% Erode main image
xdim = size(im,2); ydim = size(im,1);
no = wiener2(im,[5 5]);
level = graythresh(im);
bw = im2bw(im,level);
f = imfill(bw,'holes');

% Find connected clusters of pixels and perimeters
er = imerode(f,strel('disk',10));
er = imerode(er,strel('disk',4));


% % Iterative erode. Unnecessary, I think. Just once more and it seems to do
% % very well and counting the number of eggs for you.
% [L, numeggs] = bwlabel(er);

% i=1;
% while numeggs<numeggs
%     i = i+1;
%     er = imerode(er,strel('disk',10));
%     [L, numeggs] = bwlabel(er);
%     subplot(2,2,3); imshow(er);
%     pause(.1)
% end


[L, numeggs] = bwlabel(er);

% % Display
% close all;
% figure('units','normalized','position',[.1 .1 .8 .8]);
% subplot(2,2,1); imshow(im,[0 10])
% subplot(2,2,2); imshowt(f)

%% Rotate mask by 10deg intervals and convolute with the whole image, find
% local peaks to compare with CoMs
peaks=zeros(numeggs,3,36);
for j = 1:36
    close all;
    rots(j) = (j-1)*10;
    eggr = double(imrotate(egg,rots(j))); % 10 deg rotation increments
    %         G = fspecial('gaussian',[10 10],8);
    %         eggr = imfilter(eggr,G,'same');
    eggr = eggr - mean(mean(eggr));
    
    C = conv2(double(im),double(eggr),'same');
    % Find peaks iteratively
    rad = min(size(egg,1),size(egg,2))*.5;
    tmppeaks = findpeaksT(C,numeggs,rad);
    peaks(:,:,j) = tmppeaks;
    if mod(j,6)==0 disp(num2str(j)); end;
end

%% Find best convolution match and rotation for each egg
% Compare peaks to CoM's of each blob/egg. Points must be within 1/10 the
% radius to qualify. Then find the point with the highest score
best_matches=zeros(numeggs,4);

for i = 1:numeggs
    [y,x] = find(L==i);
    CoMs(i,1) = mean(x);
    CoMs(i,2) = mean(y);
end
CoMs = sortrows(CoMs,1); % sort left to right

for i = 1:numeggs
    clear scores
    
    % Calculate distances in a 2d matrix of dims [peaknum,rotnum]
    dists = squeeze(sqrt((peaks(:,1,:)-CoMs(i,1,:)).^2+(peaks(:,2,:)-CoMs(i,2,:)).^2));
    [peaknums,rotnums] = find(dists<rad/5);
    for j = 1:size(peaknums,1)
        scores(j) = peaks(peaknums(j),3,rotnums(j));
    end
    bestind = find(scores==max(scores));
    best_matches(i,:) = [peaks(peaknums(bestind),:,rotnums(bestind)) (rotnums(bestind)-1)*10];
    % columns: [x y score rot]
    
    % Plot best matches over the image
    mar = double(imrotate(ma,(rotnums(bestind)-1)*10));
    mapr = imrotate(map,(rotnums(bestind)-1)*10);
    [my,mx] = find(mapr);
    mx = mx + best_matches(i,1) - round(size(mapr,2)/2);
    my = my + best_matches(i,2) - round(size(mapr,2)/2);
    matchcoords{i} = [mx my];
    
    % Make masks for active contours
    ind1 = (best_matches(i,2)-round(size(mar,1)/2)+1):(best_matches(i,2)+round(size(mar,1)/2)-1);
    ind2 = (best_matches(i,1)-round(size(mar,2)/2)+1):(best_matches(i,1)+round(size(mar,2)/2)-1);
    acmasks{i} = zeros(size(im,1),size(im,2));
    acmasks{i}(ind1,ind2) = mar;
end


