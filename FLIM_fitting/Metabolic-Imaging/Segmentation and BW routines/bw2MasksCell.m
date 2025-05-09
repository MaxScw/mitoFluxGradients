function masks = bw2MasksCell(bw)
% Simple function to take a single bw mask of blobs. It sorts them into
% individual elements of a cell array, in order of x-center-of-mass

bwSnipBool = 0;
[ydim,xdim] = size(bw);

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
    
    [L2 num2] = bwlabel(snipmask);
    for b2 = 1:max(L2(:))
        ind2 = find(L2==b2);
        tmp2 = zeros(ydim,xdim);
        tmp2(ind2) = 1;
        % Populate final masks cell
        masks = [masks boolean(tmp2)];
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