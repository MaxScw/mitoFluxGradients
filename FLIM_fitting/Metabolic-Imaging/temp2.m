clear all;
IllProfCal = tiffread2('C:\Dropbox\data\s1_a1_zscan\cal_files\IllProfCal_40Xobj_1p2NA_1p28BHzoom.tif')
IllProfCal = double(IllProfCal.data);
Gblur=15;
for i = 1:31
%     i = 1
    im = tiffread2(['C:\Dropbox\data\s1_a1_zscan\sorted_sdts\IntTiffs_Pos0_s1_a1_zscan\fr' num2str(i,'%05i') '.tif']);
    im = double(im.data);
    im=im(:,513:end);
    % I found IllProfCal scaling often caused artifact. Edge background
    % gets artificially elevated and confused for cytoplasm. So just don't
    % do it.
%     im = im./IllProfCal.*mean(IllProfCal(:)); 
    subplot(2,2,1); imshow(im,[0 10]);
    bim = gaussianbpf(im,2,25);
    subplot(2,2,2); imshow(bim,[0 10]);
    level = graythresh(bim); BW = imbinarize(bim,level*.5);
    BW = bwareafilt(BW,[50 10^8]);
    subplot(2,2,3); imshow(BW);
    BW = imdilate(BW,strel('disk',10));
    [y,x] = find(bwperim(BW));
    subplot(2,2,4); imshow(im,[0 10]); hold on; plot(x,y,'b.')
end