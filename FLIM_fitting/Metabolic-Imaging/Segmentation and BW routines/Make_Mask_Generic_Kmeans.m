function Make_Mask_Generic_Kmeans(acqpath,poss,area_cuts,Gblur,nClusts,KeepClusts,thresh)
% Once redox side-by-side tiff intensity images have been made, identify
% the right regions to include in a decay by using kmeans to find
% similar-valued pixels.
% Inputs:
% -area_cuts: exclude regions larger or smaller than [low high]
% -Gblur: before kmeans, perform a gaussian blur on the image. Gblur is the
%   radius and sigma of the blur
% -Nclusts: number of clusters for kmeans
% -KeepClusts: index of clusters to include (ex, 3 clusters, [2 3] would
%   keep the brightest two clusters of pixels)
% -thresh: upper and lower threshold values for NADH and FAD [Nmin Nmax Fmin Fmax]

% clear all;
% acqpath = 'Z:\Lab\Tim\2015-08-11 Hembryo\em1\';
% nClusts = 4; KeepClusts = [2 3];
% Gblur = 10;
% area_cuts = [100 10000]

if acqpath(end)~='\' acqpath = [acqpath '\']; end;
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end
NADHBool = 1; FADBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [1000 10000]; end;
if ~exist('Gblur')|Gblur==-1 Gblur = 6; end;
if ~exist('nClusts')|nClusts==-1 nClusts = 3; end;
if ~exist('KeepClusts')|KeepClusts==-1 KeepClusts = [2 3]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('thresh')|thresh==-1 thresh = [0 10^6 0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
ProcDisp = 0;
close all;
G = fspecial('gaussian',[Gblur Gblur],Gblur);
Gt = fspecial('gaussian',[5 5],5);

% Check for IllProfCal.mat file
ProfBool = 0;
if exist([acqpath 'cal_files'])==7 % if new cal_files folder exists
    Dprof = dir([acqpath 'cal_files\IllProfCal*.mat']);
    if ~isempty(Dprof)
        load([acqpath 'cal_files\' Dprof(1).name]);
        IllProfCal = double(IllProfCal);
        IllProfCaldual = [IllProfCal IllProfCal];
        ProfBool = 1;
    end
else
    Dprof = dir([UpOneDir(acqpath) '\DailyFiles\IllProfCal.mat']);
    if ~isempty(Dprof)
        load([UpOneDir(acqpath) '\DailyFiles\IllProfCal.mat']);
        IllProfCal = double(IllProfCal);
        IllProfCaldual = [IllProfCal IllProfCal];
        ProfBool = 1;
    end
end
slashes = strfind(acqpath,'\');
Run = acqpath(slashes(end-1)+1:end-1);
Dpos = dir([acqpath 'sorted_sdts\*Pos*']);
Dpos(~[Dpos.isdir]) = [];
remove = [];
for i = 1:length(Dpos)
    if length(Dpos(i).name)>5
        remove = [remove i];
    end
end
Dpos(remove)=[];
if ~isempty(Dpos)
    srtpath = [acqpath 'sorted_sdts\'];
    for i = 1:size(Dpos,1)
        FLIMagepaths{i} = [srtpath 'FLIMageTiffs_' Dpos(i).name '_' Run '\'];
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
%         Dtifs{i} = dir([tifpaths{i} '*.tif']);
        DFLIMages{i} = dir([tifpaths{i} '*.tif']);
        % save figures with overlaid ROIs to do quick checks after batch processing
        ROIpaths{i} = [srtpath 'ROIsCheckGen_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = acqpath;
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
    ts = unique(sort([nameinds{PosInd,3}]))+1;
    Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if ~exist('rang') rang = [min(Zs) max(Zs)]; end;
    
    set(gcf,'PaperPositionMode', 'auto')
    clear Nmasks Fmasks Masks EggVals NADim FADim
    
    
    for i = 1:length(frames)
        if frames(i)==15
            0
        end
        t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),2}]))+1;
        Z(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),5}]))+1;
        disp(['Pos' uManPos ', Fr ' num2str(frames(i))])
        
        im = double(imread([tifpaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']));
        if ProfBool
            im = im./IllProfCaldual(:,1:size(im,2)).*mean(IllProfCaldual(:));
        end
%         FLim = imread([FLIMagepaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']);
        xdim = size(im,2); ydim = size(im,1);
        
        if NADHBool&FADBool
            xdim1ch = xdim/2;
            NADim = im(:,1:xdim1ch); %NADFLim = FLim(:,1:xdim1ch);
            Ngim = imfilter(NADim,Gt,'same');
            FADim = im(:,xdim1ch+1:end); %FADFLim = FLim(:,xdim1ch+1:end);
            Fgim = imfilter(FADim,Gt,'same');
            threshpix = [find(NADim<thresh(1)|NADim>thresh(2)); find(FADim<thresh(3)|FADim>thresh(4))];
            NADim(threshpix) = 0; FADim(threshpix) = 0;
            KMim = [NADim FADim];% KMFLim = [NADFLim FADFLim];
        elseif NADHBool&~FADBool
            NADim = im; %NADFLim = FLim;
            Ngim = imfilter(NADim,Gt,'same');
            FADim = [];
            threshpix = [find(NADim<thresh(1)|NADim>thresh(2))];
            KMim = [NADim]; %KMFLim = [NADFLim];
        elseif ~NADHBool&FADBool
            NADim = [];
            FADim = im;% FADFLim = FLim;
            Fgim = imfilter(FADim,Gt,'same');
            threshpix = [find(FADim<thresh(1)|FADim>thresh(2))];
            KMim = [FADim]; %KMFLim = [FADFLim];
        else
            error('Something wrong with channels')
        end
        
        [masks singlemask numeggs(i)] = Masks_Kmeans_FLIMages(KMim,KMim,Gblur,nClusts,KeepClusts,area_cuts,-1);
        singlemask(threshpix) = 0;
        StIms(frames(i)).data = uint8(singlemask.*86);
        
        % Bad frames?
        if isempty(find(singlemask))
            Masks{frames(i)}(1).NL = []; Masks{frames(i)}(1).FL = [];
            Masks{frames(i)}(1).NLper = []; Masks{frames(i)}(1).FLper = [];
            imshow(im,[min(im(:)) max([min(im(:))+1,max(im(:))*.4])],'Border','tight'); hold on;
            set(gcf,'PaperPositionMode', 'auto')
            text(xdim/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',25)
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            if ~ProcDisp
                saveas(imh,[ROIpaths{posnum} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
            end
        else
            imshow(im,[min(im(:)) max([min(im(:))+1,max(im(:))*.4])],'Border','tight'); hold on;
            set(gcf,'PaperPositionMode', 'auto')
            if NADHBool&FADBool
                [y,x] = find(bwperim(singlemask));
                plot(x,y,'b.','markersize',3);
                plot(x+xdim/2,y,'g.','markersize',3)
            elseif NADHBool&~FADBool
                [y,x] = find(bwperim(singlemask));
                plot(x,y,'b.','markersize',3);
            elseif ~NADHBool&FADBool
                [y,x] = find(bwperim(singlemask));
                plot(x,y,'g.','markersize',3);
            end
            %     set(gca,'position',[.52 .02 .45 .45]) % for multi-plot
            hold off;
            if ProcDisp
                pause(1)
            else
                set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
                saveas(imh,[ROIpaths{posnum} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
                pause(.01)
            end
            Masks{frames(i)}(1).NL = singlemask; Masks{frames(i)}(1).FL = singlemask;
            Masks{frames(i)}(1).NLper = bwperim(singlemask); Masks{frames(i)}(1).FLper = bwperim(singlemask);
            pause(.1)
        end
    end
    save([acqpath acqpath(slashes(end-1)+1:end-1) '_' Dpos(posnum).name '_SingleMasks.mat'],'frames','Masks');
    StkWrite(StIms, [srtpath 'Genmasks_' Dpos(posnum).name '_' Run '.tif']);
end

