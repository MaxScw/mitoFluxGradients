function [acs acscoords] = EllMasks_TouchingEggs_ActiveConts(im,acmasks,niter,sm,contbias,egg)
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

if ~exist('niter') niter = 100; end;
if ~exist('sm') sm = 2; end;
if ~exist('contbias')|contbias==-1 contbias = -.2; end;
if exist('egg')
    ind = egg;
else
    ind = 1:size(acmasks,2);
end

%% Get rid of some noise
no = wiener2(im,[5 5]);
G = fspecial('gaussian',[20 20],20);
gim = imfilter(double(no),G,'same');
%% Lastly, make it a tight boundary with an active contour
acmasks0 = acmasks;
for i = ind
    overlap = zeros(size(im,1),size(im,2));
    for k = 1:size(acmasks,2)
        if i~=k
            acmasks{i}(find(acmasks0{i}==acmasks0{k}))=0;
            %             acmasks{i}(find(acmasks0{i}==imerode(acmasks0{k},strel('disk',2))))=0;
        end
    end
    ind = find(overlap);
    acmasks{i}(ind)=0;
%     pause(0.5)
%     bw = activecontour(gim, imerode(acmasks{i},strel('disk',2)), 200, 'edge','SmoothFactor',1.5);
    acs{i} = activecontour(gim,imerode(acmasks{i},strel('disk',2)),niter,'Chan-Vese','SmoothFactor',sm,'ContractionBias',contbias);%Chan-Vese
    [y,x] = find(bwperim(acs{i}));
    acscoords{i} = [x y];
end

