
function [blobs bip bw L] = EllMask_Iner(im,circ_cut,egg_diam,level_fact)
% Take in an image, bpass it, threshold, find pixel clusters. Then filter
% blobs according to size and circularity, using moment of inertia tensor
% Return the mask of all circle-like structures.

xdim = size(im,2); ydim = size(im,1);
if ~exist('egg_diam') egg_diam = ydim/3; end
if ~exist('circ_cut') circ_cut = 1; end

area_cut = pi*(egg_diam/2/4)^2;
bip = bpassTS(im,ydim/50,egg_diam);
level = graythresh(bip);
bw = im2bw(bip,level*level_fact);
bw = imfill(bw,'holes');

% Find connected clusters of pixels and perimeters
[L, num] = bwlabel(bw);

% Find perimeter of pixel clusters
bwper = bwperim(bw);
blobs = [];
for i = 1:num
    % For each blob, filter for size and anisotropy, then make a cell with
    % all blobs and another with the corresponding perimeters
    temp = zeros(ydim,xdim);
    temp(find(L==i)) = 1;
    blob.L = temp;
    [blob.y,blob.x] = find(L==i);
    blob.cx = mean(blob.x);
    blob.cy = mean(blob.y);
    blob.I = [[sum((blob.y-blob.cy).^2) ...
        -sum((blob.x-blob.cx).*(blob.y-blob.cy))];...
        [-sum((blob.x-blob.cx).*(blob.y-blob.cy)) ...
        sum((blob.x-blob.cx).^2)]];
    blob.eigvals = eigs(blob.I);
    blob.eigratio = blob.eigvals(1)/blob.eigvals(2);
    blob.circularity = abs(1-blob.eigratio);
    
    % Filter:
    if blob.circularity<circ_cut & size(blob.x,1)>area_cut
        blob.Lper = bwlabel(bwperim(blob.L));
        [blob.yper,blob.xper] = find(blob.Lper);
        blobs = [blobs blob];
    end
end
disp('')






