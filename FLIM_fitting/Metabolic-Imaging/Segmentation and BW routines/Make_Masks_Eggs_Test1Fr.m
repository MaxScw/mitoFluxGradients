function Make_Masks_Eggs_Test1Fr(path,poss,area_cuts,Gblur,nClusts,KeepClusts,thresh,circ_cut,frames)
% Once redox side-by-side tiff intensity images have been made, calculate
% the redox ratio from the intensities of the eggs. Find circles using the
% bpass, thresh, then calculate principal axes of the moment of inertia
% tensor to judge circularity. Then sum intensity only within those circles.


% Versions:
% 2016-5-09: Taken from Make_Masks_Eggs_DBTrack. Run on a single frame to
% see if mask params are working as intended.

% clear all;
% path = 'Z:\Lab\Tim\2015-11-05 Ibidi Temp Cals\KCN_FCCP_Cont\';
% poss = 0;


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
if ~exist('frames')|frames==-1 frames = 1; end % Unless specified, assume roughly round eggs
load([path '\name_indexes.mat']);
ProcDisp = 0;
close all;
G = fspecial('gaussian',[Gblur Gblur],Gblur);
Gt = fspecial('gaussian',[5 5],5);

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
        FLIMagepaths{i} = [srtpath 'FLIMageTiffs_' Dpos(i).name '_' Run '\'];
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
%         Dtifs{i} = dir([tifpaths{i} '*.tif']);
        DFLIMages{i} = dir([tifpaths{i} '*.tif']);
        % save figures with overlaid ROIs to do quick checks after batch processing
        ROIpaths{i} = [srtpath 'ROIsCheck_' Dpos(i).name '_' Run '\'];
%         [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = path;
end

for posnum = 1:size(Dpos,1)
    uManPos = Dpos(posnum).name(4:end);
%     strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
%     uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    PosInd = strcmp(nameinds(:,2),uManPos); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),uManPos)
            continue;
        end
    end
    
%     frames = unique(sort([nameinds{PosInd&[nameinds{:,7}]>-1,6}]));
%     frames = frames(frames>0);
    ts = unique(sort([nameinds{PosInd,3}]))+1;
    Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if ~exist('rang') rang = [min(Zs) max(Zs)]; end;
    
    clear Nmasks Fmasks Masks Masks0 EggVals
    
    
    clear NADim FADim
    Coords = [];
    CumDriftX = 0; CumDriftY = 0; CumDrifts = [];
    for i = 1:size(frames,2)
        imh = figure;
        
        t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),3}]))+1;
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
            Ngim = imfilter(NADim,Gt,'same');
            FADim = im(:,xdim1ch+1:end)./IllProfCal(:,1:xdim1ch).*mean(IllProfCal(:)); %FADFLim = FLim(:,xdim1ch+1:end);
            Fgim = imfilter(FADim,Gt,'same');
            threshpix = [find(Ngim<thresh(1)|Ngim>thresh(2)); find(Fgim<thresh(3)|Fgim>thresh(4))];
            NADim(threshpix) = 0; FADim(threshpix) = 0;
            KMim = [NADim FADim]; %KMFLim = [NADFLim FADFLim];
        elseif NADHBool&~FADBool
            xdim1ch = xdim;
            NADim = im./IllProfCal.*mean(IllProfCal(:)); %NADFLim = FLim;
            Ngim = imfilter(NADim,Gt,'same');
            FADim = [];
            threshpix = [find(Ngim<thresh(1)|Ngim>thresh(2))];
            KMim = [NADim]; %KMFLim = [NADFLim];
        elseif ~NADHBool&FADBool
            xdim1ch = xdim;
            NADim = [];
            FADim = im./IllProfCal.*mean(IllProfCal(:)); %FADFLim = FLim;
            Fgim = imfilter(FADim,Gt,'same');
            threshpix = [find(Fgim<thresh(1)|Fgim>thresh(2))];
            KMim = [FADim]; %KMFLim = [FADFLim];
        else
            error('Something wrong with channels')
        end
        
        [masks singlemask numeggs(i)] = Masks_Kmeans_FLIMages(KMim,KMim,Gblur,nClusts,KeepClusts,area_cuts,-1);
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
                Masks0{frames(i)}(eg).DdCoM = [DeDriftX DeDriftY];
                Coords = [Coords; DeDriftX DeDriftY frames(i)];
            end
        else
            Masks0{frames(i)} = [];
        end
        
        % Display results
        xdim = size(KMim,2); ydim = size(KMim,1);
        imshow(KMim,[min(KMim(:)) max([max(max(KMim))*.2 (min(KMim(:))+1)])],'Border','tight'); hold on;
        for eg = 1:size(Masks0{frames(i)},2)
            if NADHBool&FADBool
                [y,x] = find(Masks0{frames(i)}(eg).NLper);
                plot(x,y,'b.','markersize',3)
                text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
                [y,x] = find(Masks0{frames(i)}(eg).FLper);
                plot(x+xdim1ch,y,'g.','markersize',3)
                text(mean(x)+xdim1ch,mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
            elseif NADHBool&~FADBool
                [y,x] = find(Masks0{frames(i)}(eg).NLper);
                plot(x,y,'b.','markersize',3)
                text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
            elseif ~NADHBool&FADBool
                [y,x] = find(Masks0{frames(i)}(eg).FLper);
                plot(x,y,'g.','markersize',3)
                text(mean(x),mean(y),num2str(eg),'HorizontalAlignment','center','fontsize',14,'color','r');
            else
                error('Check channels')
            end
        end
        set(gcf,'paperPositionMode','auto')
        set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
        set(gcf,'inverthardcopy','off')
    end
end
