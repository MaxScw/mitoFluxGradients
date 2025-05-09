function Make_Mask_Generic_WoW_AnyFLIMageSeries(path,area_cuts,Gblur,nClusts,KeepClusts,thresh,WoWBool)
% Adaptation of 'Make_Mask_Generic_WoW.m', which only works on redox images
% acquired with uManager channels. This version works on any arbitrary
% series of FLIM images saved to a folder.
% IMPORTANT: You must run 'IntensitymatrixTiffs.m' on the path of FLIMages
% first. This creates Tiffs in a subfolder called 'IntTiffs', which are
% used to do image segmentation (mask creation).
% Inputs:
% -area_cuts: exclude regions larger or smaller than [low high]
% -Gblur: before kmeans, perform a gaussian blur on the image. Gblur is the
%   radius and sigma of the blur
% -Nclusts: number of clusters for kmeans
% -KeepClusts: index of clusters to include (ex, 3 clusters, [2 3] would
%   keep the brightest two clusters of pixels)

% clear all;
% path = 'C:\Dropbox\data\2016-11-17ROS_Test_marta\better ones\';
% nClusts = 2; KeepClusts = [2];
% Gblur = 10;
% area_cuts = [100 10000]

if path(end)~='\' path = [path '\']; end;
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [1000 10^6]; end;
if ~exist('Gblur')|Gblur==-1 Gblur = 6; end;
if ~exist('nClusts')|nClusts==-1 nClusts = 2; end;
if ~exist('KeepClusts')|KeepClusts==-1 KeepClusts = [2]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('thresh')|thresh==-1 thresh = [0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('NumWells')|NumWells==-1 NumWells = 9; end % Unless specified, assume 9-well dishes. Enter '16' if you're using 16-well dishes
if ~exist('WoWBool')|WoWBool==-1 WoWBool = 1; end % Unless specified, WoW dishes. Enter '0' if imaging different type of sample
if NumWells==9
    WellRad = 210;
elseif NumWells==16
    WellRad = 140;
else
    WellRad = NumWells; % If not 9 or 16, assume custom radius entered
end
G = fspecial('gaussian',[5 5],5);


close all;

% Check for IllProfCal.mat file
ProfBool = 0;
if exist([path 'cal_files'])==7 % if new cal_files folder exists
    Dprof = dir([path 'cal_files\IllProfCal*.mat']);
    if ~isempty(Dprof)
        load([path 'cal_files\' Dprof(1).name]);
        IllProfCal = double(IllProfCal);
        IllProfCaldual = [IllProfCal IllProfCal];
        ProfBool = 1;
    end
else
    Dprof = dir([UpOneDir(path) '\DailyFiles\IllProfCal.mat']);
    if ~isempty(Dprof)
        load([UpOneDir(path) '\DailyFiles\IllProfCal.mat']);
        IllProfCal = double(IllProfCal);
        IllProfCaldual = [IllProfCal IllProfCal];
        ProfBool = 1;
    end
end

slashes = strfind(path,'\');
Run = path(slashes(end-1)+1:end-1);
tifpath = [path 'IntTiffs\'];
ROIpath = [path 'ROIsCheckGen\'];
[a,b] = mkdir(ROIpath);
Dtiffs = dir([tifpath '*.tif']);
L = size(Dtiffs,1);

imh = figure;
set(gcf,'PaperPositionMode', 'auto')

for i = 1:L
    
    im = double(imread([tifpath Dtiffs(i).name]));
    xdim = size(im,2); ydim = size(im,1);
    if WoWBool
        % Find peripheral well area and exclude
        Gwell = fspecial('gaussian',[10 10],10);
        [Y,X] = meshgrid(-ydim/2:(ydim/2-1),(-xdim/2:(xdim/2-1))); Z = X.^2+Y.^2;
        ind = find(sqrt(Z)<WellRad*.8); 
        im2 = im; im2(ind) = 0;
        gim = uint8(imfilter(InvertIm(im2),Gwell,'same'));
        CircMask = ~Make_Circ_Mask(xdim,ydim,0,0,WellRad);
        C = conv2(double(gim),double(CircMask),'same');
        [Cy,Cx] = find(C==max(C(:)));
    end
    
    disp([Dtiffs(i).name])
    
    if ProfBool
        im = im./IllProfCal.*mean(IllProfCal(:));
    end
    gim = imfilter(im,G,'same');
            
    threshpix = find(im<thresh(1)|im>thresh(2));
    im(threshpix) = 0;
    if WoWBool im = WoWCircCrop(im,Cx,Cy,WellRad); end
    %     NADFLim = WoWCircCrop(NADFLim,Cx,Cy,WellRad);
    KMim = [im]; %KMFLim = [NADFLim FADFLim];
    % Sclae image by IllProfCal image for more accurate pixel
    % segmentation
    KMim = KMim./IllProfCal.*mean(IllProfCal(:));
    
    [masks singlemask numeggs(i)] = Masks_Kmeans_FLIMages(KMim,-1,Gblur,nClusts,KeepClusts,area_cuts,-1);
    singlemask(threshpix) = 0;
    StIms(frames(i)).data = uint8(singlemask.*86);
    
    % Bad frames?
    if isempty(find(singlemask))
        Masks{i}(1).NL = []; Masks{i}(1).FL = [];
        Masks{i}(1).NLper = []; Masks{i}(1).FLper = [];
        imshow(im,[min(im(:)) max([min(im(:))+1,max(im(:))*.4])],'Border','tight'); hold on;
        set(gcf,'PaperPositionMode', 'auto')
        text(xdim/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',25)
        set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
        saveas(imh,[ROIpath Dtiffs(i).name '.png'],'png');
    else
        imshow(im,[min(im(:)) max([min(im(:))+1,max(im(:))*.4])],'Border','tight'); hold on;
        set(gcf,'PaperPositionMode', 'auto')
        [y,x] = find(bwperim(singlemask));
        plot(x,y,'r.','markersize',3);
        %     set(gca,'position',[.52 .02 .45 .45]) % for multi-plot
        hold off;
        pause(1)
        set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
        saveas(imh,[ROIpath Dtiffs(i).name '.png'],'png');
        pause(.01)
        Masks{i}(1).NL = singlemask; Masks{i}(1).FL = singlemask;
        Masks{i}(1).NLper = bwperim(singlemask); Masks{i}(1).FLper = bwperim(singlemask);
        pause(.1)
    end
end
save([path path(slashes(end-1)+1:end-1) '_SingleMasks.mat'],'Masks');
StkWrite(StIms, [path 'IntTiffs\Genmasks.tif']);

