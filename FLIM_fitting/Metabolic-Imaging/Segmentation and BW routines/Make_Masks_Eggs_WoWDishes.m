function Make_Masks_Eggs_WoWDishes(path,poss,area_cuts,Gblur,nClusts,KeepClusts,thresh,circ_cut,NumWells)
% Once redox side-by-side tiff intensity images have been made, calculate
% the redox ratio from the intensities of the eggs. Find circles using the
% bpass, thresh, then calculate principal axes of the moment of inertia
% tensor to judge circularity. Then sum intensity only within those circles.
% Adapted from 'Make_Masks_Eggs_DBTrack.m', but with added step of
% identifying the peripheral well fluorescence and excluding it (set to 0)
% before trying to find the embryos.

% Versions:
% 2015-10-23: Add drift compensation by cross correlating successive images
%  with fft2's. Incorporate drifts in coords before using 'track'
% 2015-09-22: CANCEL - Add rotation and drift adjustment using CoM and moment of
% inertia tensor of the whole image. Change egg tracking to do a local
% search around where you expect the egg to be.
% 2015-07-17: Added a radius filter once eggs are detected. Blocks out
% pixels outside of a certain radius of eggs. Otherwise, bright crap-chunks
% can throw off the k-means
% 2015-07-13: Switched to k-means clustering and made that
% 'Make_Masks_Eggs'. See 'Make_Masks_Eggs_Erode_ActConts' for previous
% method of 'EllMask_blur_fillholes' and active contours.
% 2015-02-09: Added NADHBool and FADBool so it could work for dual or not
% dual channel acquisitions. Name change from
% 'Make_Masks_Eggs_Iner_DualMasks' to 'Make_Masks_Eggs'
% 2014-10-17: nameinds, Mask cells, decay structs, and fit structs have
% been re-indexed to align better. This is updated to match the new scheme.
% Also add bead masks if they are there.
% -Before: Use iterative approach to find the most eggs (or use user eggs).
% 'Make_Masks_Eggs_Iner_NADMask' uses the NAD mask for both channels. This
% version fits respective channels and saves them as 'NL' and 'FL' fields

% clear all;
% path = 'C:\Users\Tim\Documents\Academic - Research\Data\Dish1_oldyoung\';
% poss = 0;

DeDriftBool = 0; 
% Wrote this dedrifting provision for when samples drifted more, but these
% days I use WoW dishes, so no collective drift occurs, so don't use.

if path(end)~='\' path = [path '\']; end;
load([path '\name_indexes.mat']);
NADHBool = 1; FADBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [1000 10000]; end;
if ~exist('Gblur')|Gblur==-1 Gblur = 6; end;
if ~exist('nClusts')|nClusts==-1 nClusts = 3; end;
if ~exist('KeepClusts')|KeepClusts==-1 KeepClusts = [2 3]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('thresh')|thresh==-1 thresh = [0 10^6 0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('circ_cut')|circ_cut==-1 circ_cut = 10; end % Unless specified, assume roughly round eggs
if ~exist('NumWells')|NumWells==-1 NumWells = 9; end % Unless specified, assume 9-well dishes. Enter '16' if you're using 16-well dishes
if NumWells==9
    WellRad = 210;
elseif NumWells==16
    WellRad = 140;
else
    WellRad = NumWells; % If not 9 or 16, assume custom radius entered
end
G = fspecial('gaussian',[5 5],5);

load([path '\name_indexes.mat']);
ProcDisp = 0;
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
Dpos = dir([path 'sorted_sdts\*Pos*']);
Dpos(~[Dpos.isdir]) = [];
remove = [];
for i = 1:length(Dpos)
    if ~strcmp(Dpos(i).name(1:3),'Pos')
        remove = [remove i];
    end
end
Dpos(remove)=[];
if ~isempty(Dpos)
    srtpath = [path 'sorted_sdts\'];
    for i = 1:size(Dpos,1)
%         FLIMagepaths{i} = [srtpath 'FLIMageTiffs_' Dpos(i).name '_' Run '\'];
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
%         Dtifs{i} = dir([tifpaths{i} '*.tif']);
%         DFLIMages{i} = dir([tifpaths{i} '*.tif']);
        % save figures with overlaid ROIs to do quick checks after batch processing
        ROIpaths{i} = [srtpath 'ROIsCheck_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = path;
end
imh = figure;

for posnum = 1:size(Dpos,1)
    uManPos = Dpos(posnum).name(4:end);
%     strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
%     uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    PosInd = strcmp(nameinds(:,3),uManPos); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),uManPos)
            continue;
        end
    end
    
    frames = unique(sort([nameinds{PosInd&[nameinds{:,7}]>-1,6}]));
    frames = frames(frames>0);
    ts = unique(sort([nameinds{PosInd,2}]))+1;
    Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if ~exist('rang') rang = [min(Zs) max(Zs)]; end;
    
    set(gcf,'PaperPositionMode', 'auto')
    clear Nmasks Fmasks Masks Masks0 EggVals
    
    
    clear NADim FADim
    Coords = [];
    CumDriftX = 0; CumDriftY = 0; CumDrifts = [];
    
    % Find peripheral well area and exclude
    im = double(imread([tifpaths{posnum} 'fr' num2str(frames(1),'%05i') '.tif']));
    xdim1ch = size(im,2); ydim = size(im,1);
    if NADHBool&FADBool
        xdim1ch = xdim1ch/2;
        im1ch = im(:,xdim1ch+1:end);
    elseif xdim1ch==ydim
        xdim1ch = xdim1ch;
        im1ch = im;
    else
        error('Something wrong with image dimensions.')
    end
    Gwell = fspecial('gaussian',[10 10],10);
    [Y,X] = meshgrid(-ydim/2:(ydim/2-1),(-xdim1ch/2:(xdim1ch/2-1))); Z = X.^2+Y.^2;
    ind = find(sqrt(Z)<WellRad*.8);
    im1ch_2 = im1ch; im1ch_2(ind) = 0;
    gim = uint8(imfilter(InvertIm(im1ch_2),Gwell,'same'));
    CircMask = ~Make_Circ_Mask(xdim1ch,ydim,0,0,WellRad);
    %         CircMask = double(CircMask)*90+10;
    C = conv2(double(gim),double(CircMask),'same');
    [Cy,Cx] = find(C==max(C(:)));
    
    for i = 1:length(frames)
        
        t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),2}]))+1;
        Z(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),5}]))+1;
        disp(['Pos' uManPos ', Fr ' num2str(frames(i))])
        if i > 1 PrevIm = im; end
        im = double(imread([tifpaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']));
        
        if i > 1
            xcr = fftshift(ifft2(fft2(PrevIm).*conj(fft2(im))));
            [DriftY DriftX] = find(xcr==max(xcr(:))); 
            DriftX = round(DriftX-xdim/2); DriftY = round(DriftY-ydim/2);
            CumDriftX = CumDriftX + DriftX(1); CumDriftY = CumDriftY  + DriftY(1); CumDrifts = [CumDrifts;[CumDriftX CumDriftY]];
        end
        
%         FLim = imread([FLIMagepaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']);
        xdim = size(im,2); ydim = size(im,1);
        %         if ~exist('feat_diam') feat_diam = xdim; end;
        if NADHBool&FADBool
            xdim1ch = xdim/2;
            % Load ims and scale image by IllProfCal image for more accurate pixel segmentation
            NADim = im(:,1:xdim1ch)./IllProfCal(:,1:xdim1ch).*mean(IllProfCal(:)); %NADFLim = FLim(:,1:xdim1ch);
            Ngim = imfilter(NADim,G,'same');
            FADim = im(:,xdim1ch+1:end)./IllProfCal(:,1:xdim1ch).*mean(IllProfCal(:)); %FADFLim = FLim(:,xdim1ch+1:end);
            Fgim = imfilter(FADim,G,'same');
            threshpix = [find(Ngim<thresh(1)|Ngim>thresh(2)); find(Fgim<thresh(3)|Fgim>thresh(4))];
            NADim(threshpix) = 0; FADim(threshpix) = 0;
            NADim = WoWCircCrop(NADim,Cx,Cy,WellRad); FADim = WoWCircCrop(FADim,Cx,Cy,WellRad);
%             NADFLim = WoWCircCrop(NADFLim,Cx,Cy,WellRad); %FADFLim = WoWCircCrop(FADFLim,Cx,Cy,WellRad);
            KMim = [NADim FADim]; %KMFLim = [NADFLim FADFLim];
            % Sclae image by IllProfCal image for more accurate pixel
            % segmentation
        elseif NADHBool&~FADBool
            xdim1ch = xdim;
            NADim = im./IllProfCal.*mean(IllProfCal(:)); %NADFLim = FLim;
            Ngim = imfilter(NADim,G,'same');
            FADim = [];
            threshpix = [find(Ngim<thresh(1)|Ngim>thresh(2))];
            NADim(threshpix) = 0;
            NADim = WoWCircCrop(NADim,Cx,Cy,WellRad); 
            KMim = [NADim]; %KMFLim = [NADFLim];
        elseif ~NADHBool&FADBool
            xdim1ch = xdim;
            NADim = [];
            FADim = im./IllProfCal.*mean(IllProfCal(:)); %FADFLim = FLim;
            Fgim = imfilter(FADim,G,'same');
            threshpix = [find(Fgim<thresh(1)|Fgim>thresh(2))];
            FADim(threshpix) = 0;
            FADim = WoWCircCrop(FADim,Cx,Cy,WellRad); 
            KMim = [FADim]; %KMFLim = [FADFLim];
        else
            error('Something wrong with channels')
        end


        [masks singlemask numeggs(i)] = Masks_Kmeans_FLIMages(KMim,-1,Gblur,nClusts,KeepClusts,area_cuts,-1);
        singlemask(threshpix) = 0;
        for mn = 1:size(masks,2) masks{mn}(threshpix) = 0; end
        
        % Get Fmasks if dual channel
        if NADHBool&FADBool
            Nmasks = masks;
            Fmasks = masks;
        elseif NADHBool&~FADBool
            Fmasks = [];
        elseif ~NADHBool&FADBool
            Fmasks = masks;
            Nmasks = [];
        else
            error('Something wrong with channels')
        end
        
        % FINAL MASKS. Plug them into the struct
        if ~isempty(masks)
            for eg = 1:size(masks,2)
                Masks0{frames(i)}(eg).NL = masks{eg}; Masks0{frames(i)}(eg).FL = masks{eg};
                Masks0{frames(i)}(eg).NLper = bwperim(masks{eg}); Masks0{frames(i)}(eg).FLper = bwperim(masks{eg});
                [y,x] = find(masks{eg}); Masks0{frames(i)}(eg).CoM = [mean(x) mean(y)];
                DeDriftX = mean(x)+CumDriftX; DeDriftY = mean(y)+CumDriftY;
                if DeDriftBool
                    Masks0{frames(i)}(eg).DdCoM = [DeDriftX DeDriftY];
                    Coords = [Coords; DeDriftX DeDriftY frames(i)];
                else
                    Masks0{frames(i)}(eg).DdCoM = [mean(x) mean(y)];
                    Coords = [Coords; mean(x) mean(y) frames(i)];
                end
            end
        else
            Masks0{frames(i)} = [];
        end
        
    end
    
    % Calculate the trajectories using Dan Blaire's simple code
    if isempty(Coords)
        continue;
    end
    param.mem = 3; param.good = 1; param.dim = 2; param.quiet = 0; 
    [trajs] = trackBlair(Coords,50,param);
    
    for i = 1:length(frames)
        im = double(imread([tifpaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']));
        % Bad frames?
        if isempty(Masks0{frames(i)})
            for eg = 1:size(Masks0{frames(i)},2)
                Masks{frames(i)}(eg).NL = []; Masks{frames(i)}(eg).FL = [];
                Masks{frames(i)}(eg).NLper = []; Masks{frames(i)}(eg).FLper = [];
            end
            imshow(im,[min(min(im)) max(max(im))*.4],'Border','tight'); hold on;
            set(gcf,'PaperPositionMode', 'auto')
            text(xdim/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',25)
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            if ~ProcDisp
                saveas(imh,[ROIpaths{posnum} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
            end
            nameinds((PosInd)&([nameinds{:,6}]==frames(i)),[7 8 9]) = num2cell(-1);
            continue;
        end
        
        % Resort the structure egg elements according to particle IDs
        for eg = 1:size(Masks0{frames(i)},2)
            trajID = trajs(trajs(:,1)==Masks0{frames(i)}(eg).DdCoM(1)&trajs(:,2)==Masks0{frames(i)}(eg).DdCoM(2)&trajs(:,3)==frames(i),4);
            Masks{frames(i)}(trajID) = Masks0{frames(i)}(eg);
        end
        
        if ProfBool
            im = im./IllProfCaldual(:,1:size(im,2)).*mean(IllProfCaldual(:));
        end
        xdim = size(im,2); ydim = size(im,1);
        imshow(im,[min(im(:)) max([max(max(im))*.2 (min(im(:))+1)])],'Border','tight'); hold on;
        for eg = 1:size(Masks{frames(i)},2)
            if NADHBool&FADBool
                [y,x] = find(Masks{frames(i)}(eg).NLper);
                plot(x,y,'b.','markersize',3)
                text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
                [y,x] = find(Masks{frames(i)}(eg).FLper);
                plot(x+xdim1ch,y,'g.','markersize',3)
                text(mean(x)+xdim1ch,mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
            elseif NADHBool&~FADBool
                [y,x] = find(Masks{frames(i)}(eg).NLper);
                plot(x,y,'b.','markersize',3)
                text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
            elseif ~NADHBool&FADBool
                [y,x] = find(Masks{frames(i)}(eg).FLper);
                plot(x,y,'g.','markersize',3)
                text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
            else
                error('Check channels')
            end
        end
        set(gcf,'paperPositionMode','auto')
        set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
        set(gcf,'inverthardcopy','off')
        saveas(imh,[ROIpaths{posnum} 'ROI_im' num2str(frames(i) ,'%03i') '.png'],'png');
        hold off;
        pause(0.1)
    end
    
    save([path path(slashes(end-1)+1:end-1) '_' Dpos(posnum).name '_Masks.mat'],'frames','Masks','Coords','trajs');
end
save([path 'name_indexes.mat'],'nameinds')
