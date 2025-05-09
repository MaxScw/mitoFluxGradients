function Make_Mask_Generic_WoW(acqpath,poss,area_cuts,Gblur,nClusts,KeepClusts,thresh,NumWells,frames)
% Once redox side-by-side tiff intensity images have been made, identify
% the right regions to include in a decay by using kmeans to find
% similar-valued pixels.
% Adaptation of 'Make_Mask_Generic_Kmeans.m' to crop out edge of WoW dishes
% Inputs:
% -area_cuts: exclude regions larger or smaller than [low high]
% -Gblur: before kmeans, perform a gaussian blur on the image. Gblur is the
%   radius and sigma of the blur
% -Nclusts: number of clusters for kmeans
% -KeepClusts: index of clusters to include (ex, 3 clusters, [2 3] would
%   keep the brightest two clusters of pixels)

% clear all;
% acqpath = 'Z:\Lab\Marta\Racowsky Collab\2017-10-11 Cumulus cells\s1_a1\';
% poss = 1;
% nClusts = 3; KeepClusts = [2 3];
% Gblur = 5;
% area_cuts = [50 10^6];
% thresh = [0 30 0 30];
% NumWells = 0;

if acqpath(end)~='\' acqpath = [acqpath '\']; end;
load([acqpath 'multiD_indices.mat']);
NADHBool = 1; FADBool = 1; UserChBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'UserChan'))) UserChBool = 0; end
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [0 10^6]; end;
if ~exist('Gblur')|Gblur==-1 Gblur = 5; end;
if ~exist('nClusts')|nClusts==-1 nClusts = 3; end;
if ~exist('KeepClusts')|KeepClusts==-1 KeepClusts = [2 3]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('thresh')|thresh==-1 thresh = [0 10^6 0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('NumWells')|NumWells==-1 NumWells = 0; end % Unless specified, assume no well dishe. Enter '16' if you're using 16-well dishes
if NumWells==9
    WellRad = 210;
elseif NumWells==16
    WellRad = 140;
elseif NumWells==0 % If NumWells is 0, assume no wells and no cropping necessary
    WellRad = -1;
else
    WellRad = NumWells; % If not 9 or 16, assume custom radius entered
end
G = fspecial('gaussian',[Gblur Gblur],Gblur);
Gt = fspecial('gaussian',[Gblur Gblur],Gblur);

% Determine which channels are present
uchans = unique(nameinds(:,4));
if ~isempty(find(strcmp(uchans,'NADH'))) Chind(1) = 1; end
if ~isempty(find(strcmp(uchans,'FAD'))) Chind(2) = 2; end
if ~isempty(find(strcmp(uchans,'UserChan'))) Chind(3) = 3; end
Chind(Chind==0)=[];
ChLabs = {'NADH','FAD','UserChan'};

ProcDisp = 0;
close all;

% Find cal_files folder
if exist([acqpath 'cal_files'])==7 % if cal_files is in acqpath
    calpath = [acqpath 'cal_files\'];
elseif exist([UpOneDir(acqpath) 'cal_files'])==7 % if cal_files in daypath
    calpath  = [UpOneDir(acqpath) 'cal_files\'];
else
    error('Cannot locate cal_files folder. Place in daypath or acqpath, please');
end

% Load IllProfCal.mat file
ProfBool = 0;
Dprof = dir([calpath 'IllProfCal*.mat']);
if ~isempty(Dprof)
    load([calpath Dprof(1).name]);
    IllProfCal = double(IllProfCal);
    %         IllProfCaldual = [IllProfCal IllProfCal];
    ProfBool = 1;
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
%         FLIMagepaths{i} = [srtpath 'FLIMageTiffs_' Dpos(i).name '_' Run '\'];
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
        %         Dtifs{i} = dir([tifpaths{i} '*.tif']);
%         DFLIMages{i} = dir([tifpaths{i} '*.tif']);
        % save figures with overlaid ROIs to do quick checks after batch processing
        ROIpaths{i} = [srtpath 'ROIsCheck_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = acqpath;
end
imh = figure;

for pos = 1:size(Dpos,1)
    PosNum = Dpos(pos).name(4:end);
    %     strnums = sscanf(PosNum ,'%g'); %Find the numbers in the name
    %     PosNum = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),PosNum)
            continue;
        end
    end
    
    frames = unique(sort([nameinds{PosInd&[nameinds{:,7}]>-1,6}]));
    frames = frames(frames>0);
    ts = unique(sort([nameinds{PosInd,2}]))+1;
    Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if ~exist('rang') rang = [min(Zs) max(Zs)]; end
    
    set(gcf,'PaperPositionMode', 'auto')
    clear Nmasks Fmasks Umasks Masks EggVals NADim FADim
    
    % Find peripheral well area and exclude
    im = double(imread([tifpaths{pos} 'fr' num2str(frames(1),'%05i') '.tif']));
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
    % If WellRad = 0, this signifies that we aren't using WoW dishes
    if WellRad>0
        Gwell = fspecial('gaussian',[10 10],10);
        [Y,X] = meshgrid(-ydim/2:(ydim/2-1),(-xdim1ch/2:(xdim1ch/2-1))); Z = X.^2+Y.^2;
        ind = find(sqrt(Z)<WellRad*.8); 
        im1ch_2 = im1ch; im1ch_2(ind) = 0;
        gim = uint8(imfilter(InvertIm(im1ch_2),Gwell,'same'));
        CircMask = ~Make_Circ_Mask(xdim1ch,ydim,0,0,WellRad);
        %         CircMask = double(CircMask)*90+10;
        C = conv2(double(gim),double(CircMask),'same');
        [Cy,Cx] = find(C==max(C(:)));
    end
    
    for i = 1:length(frames)
        
        crpind = [];
        %         t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),2}]))+1;
        %         Z(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),5}]))+1;
        disp(['Pos' PosNum ', Fr ' num2str(frames(i))])
        
        im = double(imread([tifpaths{pos} 'fr' num2str(frames(i),'%05i') '.tif']));
        
        % Obsolete - I don't really use FLIMages anymore.
        %         FLim = imread([FLIMagepaths{pos} 'fr' num2str(frames(i),'%05i') '.tif']);
        %         FLim = imread([FLIMagepaths{pos} 'fr' num2str(frames(i),'%05i') '.tif']);
        xdim = size(im,2); ydim = size(im,1);
        
        if NADHBool&FADBool
            % Load ims and scale image by IllProfCal image for more accurate pixel segmentation
            NADim = im(:,1:xdim1ch)./IllProfCal.*mean(IllProfCal(:)); %NADFLim = FLim(:,1:xdim1ch);
            Ngim = imfilter(NADim,Gt,'same');
            FADim = im(:,xdim1ch+1:end)./IllProfCal.*mean(IllProfCal(:)); %FADFLim = FLim(:,xdim1ch+1:end);
            Fgim = imfilter(FADim,Gt,'same');
            threshpix = [find(Ngim<thresh(1)|Ngim>thresh(2)); find(Fgim<thresh(3)|Fgim>thresh(4))];
            NADim(threshpix) = 0; FADim(threshpix) = 0;
            if WellRad> 0 [NADim,crpind] = WoWCircCrop(NADim,Cx,Cy,WellRad); FADim = WoWCircCrop(FADim,Cx,Cy,WellRad); end
            %NADFLim = WoWCircCrop(NADFLim,Cx,Cy,WellRad); FADFLim = WoWCircCrop(FADFLim,Cx,Cy,WellRad);
            KMim = [NADim FADim]; %KMFLim = [NADFLim FADFLim];
        elseif NADHBool&~FADBool
            NADim = im./IllProfCal.*mean(IllProfCal(:)); %NADFLim = FLim;
            Ngim = imfilter(NADim,Gt,'same');
            FADim = [];
            threshpix = [find(Ngim<thresh(1)|Ngim>thresh(2))];
            NADim(threshpix) = 0;
            if WellRad> 0 [NADim,crpind] = WoWCircCrop(NADim,Cx,Cy,WellRad); end
            KMim = [NADim]; %KMFLim = [NADFLim];
        elseif ~NADHBool&FADBool
            NADim = [];
            FADim = im./IllProfCal.*mean(IllProfCal(:)); %FADFLim = FLim;
            Fgim = imfilter(FADim,Gt,'same');
            threshpix = [find(Fgim<thresh(1)|Fgim>thresh(2))];
            FADim(threshpix) = 0;
            if WellRad> 0 [FADim,crpind] = WoWCircCrop(FADim,Cx,Cy,WellRad); end
            KMim = [FADim]; %KMFLim = [FADFLim];
        elseif UserChBool
%             xdim1ch = xdim;
            UserChim = im./IllProfCal.*mean(IllProfCal(:)); %FADFLim = FLim;
            Ugim = imfilter(UserChim,Gt,'same');
            threshpix = [find(Ugim<thresh(1)|Ugim>thresh(2))];
            UserChim(threshpix) = 0;
            if WellRad> 0 [UserChim,crpind] = WoWCircCrop(UserChim,Cx,Cy,WellRad); end
            KMim = [UserChim]; 
        else
            error('Something wrong with channels')
        end
        
        [masks singlemask numeggs(i)] = Masks_Kmeans_FLIMages(KMim,-1,Gblur,nClusts,KeepClusts,area_cuts,-1,crpind);
        singlemask(threshpix) = 0;
        StIms(frames(i)).data = uint8(singlemask.*86);
        
        % Bad frames?
        if isempty(find(singlemask))
            Masks{frames(i)}(1).L = []; 
            Masks{frames(i)}(1).Lper = []; 
            imshow(im,[min(im(:)) max([min(im(:))+1,max(im(:))*.4])],'Border','tight'); hold on;
            set(gcf,'PaperPositionMode', 'auto')
            text(xdim/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',25)
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            if ~ProcDisp
                saveas(imh,[ROIpaths{pos} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
            end
        else
            imshow(im,[min(im(:)) max([min(im(:))+1,max(im(:))*.4])],'Border','tight'); hold on;
            set(gcf,'PaperPositionMode', 'auto')
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
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
            elseif UserChBool
                [y,x] = find(bwperim(singlemask));
                plot(x,y,'r.','markersize',3);
            else
                error('Channel problem.')
            end
            saveas(imh,[ROIpaths{pos} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
            pause(.01)
        end
        Masks{frames(i)}(1).L = singlemask; 
        Masks{frames(i)}(1).Lper = bwperim(singlemask);
        pause(.1)
    end
    save([acqpath 'JointMasks_' Dpos(pos).name '.mat'],'frames','Masks');
    StkWrite(StIms, [srtpath 'Genmasks_' Dpos(pos).name '_' Run '.tif']);

end

