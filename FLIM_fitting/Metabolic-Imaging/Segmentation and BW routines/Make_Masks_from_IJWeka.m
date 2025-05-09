function Make_Masks_from_IJWeka(path,poss,area_cuts,thresh,Gblur,SingMaskBool,MaskLab,RemDeadCells)
% Adapted from Make_Mask labeling functions (e.g. 'Make_Masks_Eggs_WoWDishes.m')
% In this version, image segmentation is performed using ImageJ (Fiji),
% using the 'Trainable Weka Segmentation v3.2.2' plugin. They are saved in
% 'sortded_sdts' as 'IJmasks_...tif' stacks. They are composite stacks of
% the masks overlaid on top of images, which are the product of NADH and
% FAD.

% Inputs:
% -area_cuts: exclude regions larger or smaller than [low high]
% -thresh: intensity thresholds - [NADHmin NADHmax FADmin FADmax]
% -SingMaskBool: 1 if you don't want individual masks, but single mask only
% -Gblur: gaussian blur
% -SingMaskBool: '1' to group all blobs into a single mask
% -MaskLab: Mask label (string)
% -RemDeadCells: '1' to remove dead regions of cells (very high in FAD and low in NADH)

% Versions:
% 2018-01-18: Changed path for IllProfCal
% 2017-04-04: Added option to do gaussian blur before threshold for
% additional selection. (Motivation was to filter out dead cells in blasts)

% clear all;
% path = 'Z:\Lab\Tim\2018-01-16 New SHG detector\s1_a1\';
% thresh = [0 20 0 20];
% SeparateMasksBool = 0;
% poss=3;
% Gblur = 5;
% RemDeadCells=1;

if path(end)~='\' path = [path '\']; end;
try     load([path 'multiD_indices.mat']); catch     load([path 'name_indexes.mat']); end
NADHBool = 1; FADBool = 1; UserChBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'UserChan'))) UserChBool = 0; end
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [0 10^6]; end;
if ~exist('Gblur')|Gblur==-1 Gblur = 1; end;
if ~exist('thresh')|thresh==-1 thresh = [0 10^6 0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('SingMaskBool')|SingMaskBool==-1 SingMaskBool = 0; end;
if ~exist('RemDeadCells')|RemDeadCells==-1 RemDeadCells = 0; end;
G = fspecial('gaussian',[Gblur Gblur],Gblur);

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
    if length(Dpos(i).name)>5
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
    clear Nmasks Fmasks Masks Masks0 EggVals NADim FADim IJfrs
    Coords = [];
    
    % Load Weka masks, label individual blobs, and track their
    % trajectoriese.
    IJmasks = tiffread2([srtpath 'IJmasks_' Dpos(posnum).name '_' Run '.tif']);
    IJmasks(2:2:end)=[];
    
    disp(['Pos' uManPos ' analyzing...'])
    
    % Note, the frames in the IJ stack correspond to the intensity tiff
    % frame numbers. Actually have to us the dir to get the right fr
    Dtiffs = dir([tifpaths{posnum} '*.tif']);
    for i = 1:length(Dtiffs) IJfrs(i) = str2num(Dtiffs(i).name(3:7)); end
    
    for i = 1:length(frames)
        
        t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),2}]))+1;
        Z(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),5}]))+1;
        
        IJfr = find(IJfrs==frames(i));
        
        im = double(imread([tifpaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']));
        %
        xdim = size(im,2); ydim = size(im,1);
        %
        if NADHBool&FADBool
            xdim1ch = xdim/2;
            NADim = im(:,1:xdim1ch)./IllProfCal.*mean(IllProfCal(:));
            Ngim = imfilter(NADim,G,'same');
            FADim = im(:,xdim1ch+1:end)./IllProfCal.*mean(IllProfCal(:));
            Fgim = imfilter(FADim,G,'same');
            AboveThresh = [find(Ngim>thresh(2)); find(Fgim>thresh(4))];
            BelowThresh = [find(Ngim<thresh(1)); find(Fgim<thresh(3))];
            threshpix = [AboveThresh;BelowThresh];
            %             threshpix = [find(NADim<thresh(1)|NADim>thresh(2)); find(FADim<thresh(3)|FADim>thresh(4))];
            NADim(threshpix) = 0; FADim(threshpix) = 0;
            %             % Scale image by IllProfCal image for more accurate pixel
            %             % segmentation
            if RemDeadCells
                [DeadMask DeadInd] = FindDeadCells(NADim,FADim);
                NADim(DeadInd) = 0; FADim(DeadInd) = 0;
            end
            ims{i} = [NADim FADim];
%         elseif NADHBool&~FADBool
%             xdim1ch = xdim;
%             NADim = im;
%             Ngim = imfilter(NADim,G,'same');
%             FADim = [];
%             AboveThresh = find(Ngim>thresh(2));
%             BelowThresh = find(Ngim<thresh(1));
%             threshpix = [AboveThresh;BelowThresh];
%             NADim(threshpix) = 0;
%             NADim = NADim./IllProfCal.*mean(IllProfCal(:));
%             ims{i} = [NADim];
%         elseif ~NADHBool&FADBool
%             xdim1ch = xdim;
%             NADim = [];
%             FADim = im;
%             Fgim = imfilter(FADim,G,'same');
%             AboveThresh = find(Fgim>thresh(2));
%             BelowThresh = find(Fgim<thresh(1));
%             threshpix = [AboveThresh;BelowThresh];
%             FADim(threshpix) = 0;
%             FADim = FADim./IllProfCal.*mean(IllProfCal(:));
%             ims{i} = [FADim];
%         elseif UserChBool
        else
            xdim1ch = xdim;
            Im = im./IllProfCal.*mean(IllProfCal(:));
            gim = imfilter(Im,G,'same');
            AboveThresh = find(gim>thresh(2));
            BelowThresh = find(gim<thresh(1));
            threshpix = [AboveThresh;BelowThresh];
            Im(threshpix) = 0;
            ims{i} = [Im];
%         else
%             error('Something wrong with channels')
        end
                
        % Build mask structure from IJmasks
        MaskVals = unique([IJmasks(IJfr).data]);
        cyto = zeros(size(IJmasks(IJfr).data,1),size(IJmasks(IJfr).data,2));
        cyto(find(IJmasks(IJfr).data==MaskVals(2)))=1; % cytoplasm
%         cyto(find(IJmasks(IJfr).data==MaskVals(3)))=1; % Bright spots
        % Use the extra 'nucleus info to exclude nuclear regions. In
        % practice, the Weka tends to under-fit the nuclear regions, so run
        % a 1-pixel dilation to make it a little bigger. Tested by eye.
        if size(MaskVals,1)>2
            NucIm = zeros(size(IJmasks(IJfr).data,1),size(IJmasks(IJfr).data,2));
            % Sometimes we have an extra class called 'brights'. Human
            % embryos have small spots with very high fluorescence. In this
            % case, the Weka trainer has (in the following order):
            % - background, cytoplasm, nucleus, brights
            % But the ordering of the pixel values aren't in increasing
            % order. Go figure. They come out:
            % - background(0), cytoplasm(86), nucleus(255), brights(125)
            % Must determine nucleus index from number of classes:
            if length(MaskVals)>3 %Brights present
                NucNum = 4;
            else
                NucNum = 3; brightind = [];
            end
            NucIm(find(IJmasks(IJfr).data==MaskVals(NucNum)))=1;
            nucind = find(imdilate(NucIm,strel('disk',1)));
            cyto(nucind)=0; 
        else
            nucind=[];
        end
        % Area filter to remove small blobs
        cyto(threshpix) = 0;
        if RemDeadCells cyto(DeadInd) = 0; end
        cyto = bwareafilt(boolean(cyto),area_cuts);
        
        singlemasks{i} = cyto;
        
%         RuthEdgeBool = 0;
%         if RuthEdgeBool 
%             FilledMask = imfill(uint8(singlemasks{i}));
%             OuterErode = imerode(FilledMask ,strel('disk',5));
%             OuterErodePix = find(FilledMask&~OuterErode);
%             singlemasks{i}(OuterErodePix) = 0;
%         end
        
        % Sort into labeled blobs
        masks = bw2MasksCell(singlemasks{i});
        
        % FINAL MASKS. Plug them into the struct
        if ~isempty(masks)
            for eg = 1:size(masks,2)
                Masks0{frames(i)}(eg).NL = masks{eg}; Masks0{frames(i)}(eg).FL = masks{eg};
                Masks0{frames(i)}(eg).NLper = bwperim(masks{eg}); Masks0{frames(i)}(eg).FLper = bwperim(masks{eg});
                [y,x] = find(masks{eg}); Masks0{frames(i)}(eg).CoM = [mean(x) mean(y)];
                Masks0{frames(i)}(eg).DdCoM = [mean(x) mean(y)];
                Coords = [Coords; mean(x) mean(y) frames(i)];
                % Looks like very bright pixels, especially in FAD channel,
                % might be indicative of apoptotic activity. So store where
                % the bright regions are for now. Maybe later can do an
                % analysis. Nuclear regions might also be useful for
                % something. Note: same array for all cells for a given
                % frame. No other obvious way to store them.
                Masks0{frames(i)}(eg).AboveThresh = AboveThresh;
                Masks0{frames(i)}(eg).BelowThresh = BelowThresh;
                Masks0{frames(i)}(eg).NucInd = nucind;
            end
        else
            Masks0{frames(i)} = [];
        end
        1;
    end
    
    % Calculate the trajectories using Dan Blaire's simple code
    if isempty(Coords)
        continue;
    end
    param.mem = 3; param.good = 1; param.dim = 2; param.quiet = 0;
    
    if length(frames)>1 
        MultiFrameBool = 1;
    else
        MultiFrameBool = 0;
    end
    if ~SingMaskBool&MultiFrameBool [trajs] = trackBlair(Coords,50,param); end
    
    for i = 1:length(frames)
        %         im = double(imread([tifpaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']));
        % Bad frames?
        if isempty(Masks0{frames(i)})
            for eg = 1:size(Masks0{frames(i)},2)
                Masks{frames(i)}(eg).NL = []; Masks{frames(i)}(eg).FL = [];
                Masks{frames(i)}(eg).NLper = []; Masks{frames(i)}(eg).FLper = [];
            end
            imshow(ims{i},[min(min(ims{i})) max(max(ims{i}))*.4],'Border','tight'); hold on;
            set(gcf,'PaperPositionMode', 'auto')
            text(xdim/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',25)
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            if ~ProcDisp
                saveas(imh,[ROIpaths{posnum} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
            end
            nameinds((PosInd)&([nameinds{:,6}]==frames(i)),[7 8 9]) = num2cell(-1);
            continue;
        end
        
        % Resort the structure egg elements according to particle IDs. Make
        % plots.
        
        for eg = 1:size(Masks0{frames(i)},2)
            if ~SingMaskBool&MultiFrameBool
                trajID = trajs(trajs(:,1)==Masks0{frames(i)}(eg).DdCoM(1)&trajs(:,2)==Masks0{frames(i)}(eg).DdCoM(2)&trajs(:,3)==frames(i),4);
                Masks{frames(i)}(trajID) = Masks0{frames(i)}(eg);
            else
                Masks = Masks0;
                trajs = -1;
            end
        end
        
        xdim = size(ims{i},2); ydim = size(ims{i},1);
        imshow(ims{i},[min(ims{i}(:)) max([max(max(ims{i}))*.2 (min(ims{i}(:))+1)])],'Border','tight'); hold on;
        for eg = 1:size(Masks{frames(i)},2)
            if NADHBool&FADBool
                [y,x] = find(Masks{frames(i)}(eg).NLper);
                plot(x,y,'b.','markersize',3)
                if ~SingMaskBool text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r'); end
                [y,x] = find(Masks{frames(i)}(eg).FLper);
                plot(x+xdim1ch,y,'g.','markersize',3)
                if ~SingMaskBool text(mean(x)+xdim1ch,mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r'); end
            elseif NADHBool&~FADBool
                [y,x] = find(Masks{frames(i)}(eg).NLper);
                plot(x,y,'b.','markersize',3)
                if ~SingMaskBool text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r'); end
            elseif ~NADHBool&FADBool
                [y,x] = find(Masks{frames(i)}(eg).FLper);
                plot(x,y,'g.','markersize',3)
                if ~SingMaskBool text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r'); end
            elseif UserChBool
                [y,x] = find(Masks{frames(i)}(eg).FLper);
                plot(x,y,'r.','markersize',3)
                if ~SingMaskBool text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r'); end
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
        
        % Make a cell of single masks as well in case you want to
        % analyze all embryos together.
        SingMasks{frames(i)}(1).NL = singlemasks{i}; SingMasks{frames(i)}(1).FL = singlemasks{i};
        SingMasks{frames(i)}(1).NLper = bwperim(singlemasks{i}); SingMasks{frames(i)}(1).FLper = bwperim(singlemasks{i});
        SingMasks{frames(i)}(1).AboveThres = Masks0{frames(i)}(1).AboveThresh;
        SingMasks{frames(i)}(1).BelowThresh = Masks0{frames(i)}(1).BelowThresh;
        SingMasks{frames(i)}(1).NucInd = Masks0{frames(i)}(1).NucInd;
    end
    
    if ~SingMaskBool save([path path(slashes(end-1)+1:end-1) '_' Dpos(posnum).name '_Masks.mat'],'frames','Masks','Coords','trajs'); end
    Masks = SingMasks;
    if exist('MaskLab')
        save([path path(slashes(end-1)+1:end-1) '_' Dpos(posnum).name '_' MaskLab '.mat'],'frames','Masks','Coords','trajs');
    else
        save([path path(slashes(end-1)+1:end-1) '_' Dpos(posnum).name '_SingleMasks.mat'],'frames','Masks','Coords','trajs');
    end
end

save([path 'name_indexes.mat'],'nameinds')
