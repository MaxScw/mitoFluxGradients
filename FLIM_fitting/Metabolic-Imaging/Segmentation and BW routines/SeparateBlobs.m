




im = NADim + FADim;
% close all
figure('position',[600 30 900 900])
subplot(2,2,1)
imshow(im,[])
imtp = imtophat(im, strel('disk', 10));
subplot(2,2,2)
imshow(imtp,[])

K = wiener2(imtp,[5 5]);
level = graythresh(K);
bw1 = im2bw(K,level*level_fact);
G = fspecial('gaussian',[15 15],10);
gim = imfilter(double(bw1),G,'same');
level = graythresh(gim);
bw = im2bw(gim,level*level_fact);
bw = imfill(bw,'holes'); 
subplot(2,2,3)
imshow(bw,[])


%%
close all
% figure('position',[600 30 900 900])
% subplot(2,2,1) 
figure;imshow(NADim + FADim,[])
D = bwdist(~bwe,'chessboard');
D = -D;
% subplot(2,2,2) 
% imshow(D,[])
D(~bwe) = -Inf;
L = watershed(D);
L(L==1)=0;
% subplot(2,2,3) 
figure;imshow(L,[])
dividers = imerode(bwe,strel('disk',2))&~L;
% subplot(2,2,4) 
figure;imshow(dividers)
% for i = 1:length(L)
%     for j = i:length(L)
%         clus1 = L(L==i);
%         clus2 = L(L==j);
%         neigh{i
        