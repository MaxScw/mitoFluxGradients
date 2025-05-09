function Make_Masks_from_IlastikProb(acqpath,poss,area_cuts,thresh,Gblur,JntMaskBool,MaskLab,MergeBool,RemDeadCells,WholeFrBool,ManROILab)
% Adapted from Make_Mask labeling functions (e.g. 'Make_Masks_Eggs_WoWDishes.m')
% In this version, image smskmentation is performed using ImageJ (Fiji),
% using the 'Trainable Weka Smskmentation v3.2.2' plugin. They are saved in
% 'sortded_sdts' as 'IJmasks_...tif' stacks. They are composite stacks of
% the masks overlaid on top of images, which are the product of NADH and
% FAD.

% Inputs:
% -area_cuts: exclude rmskions larger or smaller than [low high]
% -thresh: intensity thresholds - [NADHmin NADHmax FADmin FADmax]
% -JntMaskBool: 1 if you don't want individual masks, but single mask only
% -Gblur: gaussian blur
% -JntMaskBool: '1' to group all blobs into a single mask
% -MaskLab: Mask label (string)
% -RemDeadCells: '1' to remove dead rmskions of cells (very high in FAD and low in NADH)
% -MergeBool: '1' to look for a spreadsheet in the sample folder,
%   indicating which mask fragments should be merged (useful for blasts)
% -WholeFrBool: sometimes you just want the whole frame, e.g. taking
%   intensities of free dye solutions. Enter 1 to bypass Ilastik and take
%   the whole frame... mostly done, but not quite. Revisit when necessary.
% -ManROILab: Optional label for manual ROIs .mat file, created with 'ROIsByFrameGUI'

% % TEST in script mode:
% clear all;
% acqpath = 'C:\Dropbox\data\s1_a2_Zscans';
% ManROILab = 'ICM';
% JntMaskBool = 1;
% % thresh = [0 30 0 40];
% % poss=3;
% % Gblur = 5;
% % RemDeadCells=0;
% % MergeBool=0;
% % area_cuts=[3000 13000];

if acqpath(end)~='\' acqpath = [acqpath '\']; end;
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end
NADHBool = 1; FADBool = 1; UserChBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'UserChan'))) UserChBool = 0; end
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [00 10^6]; end;
if ~exist('Gblur')|Gblur==-1 Gblur = 1; end;
if ~exist('thresh')|thresh==-1 thresh = [0 10^6 0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
if ~exist('JntMaskBool')|JntMaskBool==-1 JntMaskBool = 1; end;
if ~exist('RemDeadCells')|RemDeadCells==-1 RemDeadCells = 0; end;
if ~exist('MergeBool')|MergeBool==-1 MergeBool = 0; end;
if ~exist('WholeFrBool')|MergeBool==-1 WholeFrBool = 0; end;
if ~exist('MaskLab')|MaskLab==-1 MaskLab = 'Masks'; end;
if ~exist('ManROILab')|ManROILab==-1 ManROILab=''; end;
% Replace MaskLab with ManROILab if it exists
MaskLab = [ManROILab 'Masks'];
G = fspecial('gaussian',[Gblur Gblur],Gblur);

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

% Check for mask_merges.xls file, used to merge blast smskments that are
% actually part of the same embryo.
if (exist([acqpath 'mask_merges.xls'])==2)&MergeBool
    [num,txt,res]= xlsread([acqpath 'mask_merges.xls']);
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
        Probpaths{i} = [srtpath 'ProbMaps_' Dpos(i).name '_' Run '\'];
        ROIpaths{i} = [srtpath 'ROIsCheck_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = acqpath;
end

for pos = 1:size(Dpos,1)
    
    PosNum = Dpos(pos).name(4:end);
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),PosNum)
            continue;
        end
    end
    
    % Check for mask merges
    if (exist([acqpath 'mask_merges.xls'])==2)&MergeBool
        for col = 1:size(res,2)
            if strcmp(num2str(res{1,col}),PosNum) ColInd = col; end
        end
        tmp = res(2:end,ColInd);
        for row = 1:length(tmp)
            if ~isnan(tmp{row})
                MergeArr{row} = str2num(tmp{row});
                MergeFrs(row) = MergeArr{row}(1);
                MergeArr{row} = MergeArr{row}(2:end);
            end
        end
        if exist('MaskLab')
            PrevMasks = load([acqpath acqpath(slashes(end-1)+1:end-1) '_' Dpos(pos).name '_' MaskLab '.mat']);
        else
            PrevMasks = load([acqpath acqpath(slashes(end-1)+1:end-1) '_' Dpos(pos).name '_Masks.mat']);
        end
        PrevMasks = PrevMasks.Masks;
    end
    
    
    frames = unique(sort([nameinds{PosInd&[nameinds{:,7}]>-1,6}]));
    frames = frames(frames>0);
    ts = unique(sort([nameinds{PosInd,2}]))+1;
    Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if ~exist('rang') rang = [min(Zs) max(Zs)]; end;
    
    clear Masks Masks0 EggVals NADim FADim IJfrs
    Coords = cell(3,1);
    
    
    disp(['Pos' PosNum ' analyzing...'])
    
    % Note, the frames in the IJ stack correspond to the intensity tiff
    % frame numbers. Actually have to us the dir to get the right fr
    Dtiffs = dir([tifpaths{pos} '*.tif']);
    for i = 1:length(Dtiffs) IJfrs(i) = str2num(Dtiffs(i).name(3:7)); end
    
    jointmasks = cell(length(frames),2);
    
    % Look for an optional polygon crop, which would be saved by
    % 'Make_Polygon_Crop.m'
    iminf = imfinfo([tifpaths{pos} 'fr' num2str(frames(1),'%05i') '.tif']);
    if exist([acqpath 'PolyCrop_Pos' num2str(pos-1) '.mat'])==2
        load([acqpath 'PolyCrop_Pos' num2str(pos-1) '.mat']);
    else
        polymask = ones(iminf.Height,iminf.Height);
    end
    crpind = find(~polymask);
    
    % Define a BG mask var, which will be the rmskions that have no
    % formskround masks in any frames.
    bgmasks{1} = ones(iminf.Height,iminf.Height); bgmasks{2} = ones(iminf.Height,iminf.Height); bgmasks{3} = ones(iminf.Height,iminf.Height);
    
    %     if ~WholeFrBool
    
    for i = 1:length(frames)
        
        %         t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),2}]))+1;
        %         Z(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),5}]))+1;
        
        if i ==17
            1;
        end
        
        % Load intensity images to plot over
        im = double(imread([tifpaths{pos} 'fr' num2str(frames(i),'%05i') '.tif']));
        xdim = size(im,2); ydim = size(im,1);
        
        % Look for an optional Manual ROI (just polygon crops for each frame
        % using 'ROIsByFrameGUI'.
        iminf = imfinfo([tifpaths{pos} 'fr' num2str(frames(1),'%05i') '.tif']);
        if exist([acqpath '\ManROIs_' ManROILab '.mat'])==2
            load([acqpath '\ManROIs_' ManROILab '.mat']);
            % Make ManROIs match this function's 'frames' array
            frs = frs(:,str2num(PosNum)+1);
            ManROIs = ManROIs(ismember(frs,frames),str2num(PosNum)+1,:);
            frs = frs(ismember(frs,frames));
            if ~isempty(ManROIs{i,1,1})
                ManROI{1} = ManROIs{i,1,1};
            else
                % Assume that in ManROI file is present, but particular
                % frame is missing, that no ROI was desired (e.g. ICM was
                % not in that frame)
                ManROI{1} = zeros(ydim,ydim);
            end
            % FAD ManROI
            if ~isempty(ManROIs{i,1,2})
                ManROI{2} = ManROIs{i,1,2};
            else
                % Assume that in ManROI file is present, but particular
                % frame is missing, that no ROI was desired (e.g. ICM was
                % not in that frame
                ManROI{2} = zeros(ydim,ydim);
            end
        else
            ManROI{1} = ones(ydim,ydim); ManROI{2} = ones(ydim,ydim);
        end
        
        % Define a new crop index, a combination of any 1st order global
        % polygon crop and any frame-by-frame polygon crops saved in the
        % folder.
        CRPind{1} = [crpind; find(~ManROI{1})];
        CRPind{2} = [crpind; find(~ManROI{2})];
        
        if NADHBool&FADBool
            xdim1ch = xdim/2;
            NADim = im(:,1:xdim1ch)./IllProfCal.*mean(IllProfCal(:));
            Ngim = imfilter(NADim,G,'same');
            FADim = im(:,xdim1ch+1:end)./IllProfCal.*mean(IllProfCal(:));
            Fgim = imfilter(FADim,G,'same');
            AboveThresh = [find(Ngim>thresh(2)); find(Fgim>thresh(4))];
            BelowThresh = [find(Ngim<thresh(1)); find(Fgim<thresh(3))];
            threshpix = [AboveThresh;BelowThresh];
            NADim(threshpix) = 0; FADim(threshpix) = 0;
            if RemDeadCells
                [DeadMask DeadInd] = FindDeadCells(NADim,FADim);
                NADim(DeadInd) = 0; FADim(DeadInd) = 0;
            end
            NADim(CRPind{1})=0; FADim(CRPind{2})=0;
            ims{frames(i)} = [NADim FADim];
        else
            xdim1ch = xdim;
            Im = im./IllProfCal.*mean(IllProfCal(:));
            gim = imfilter(Im,G,'same');
            AboveThresh = find(gim>thresh(2));
            BelowThresh = find(gim<thresh(1));
            threshpix = [AboveThresh;BelowThresh];
            Im(threshpix) = 0;
            Im(CRPind{1}) = 0;
            ims{frames(i)} = [Im];
        end
        
        % Loop over channels. 1=NADH, 2=FAD
        for ch = Chind
            % Look for NADH and FAD files in the ProbMap folder. If we find
            % either, load it into masks. Typically, just NADH will be enough,
            % but it may be necessary to sometimes use separate masks (e.g.
            % embryo drifted between NADH and FAD image)
            % NOTE: If prob maps have 4 channels, they are: mito, cyto, bg, nuc
            %       If they have 3, it is cyto, bg, nuc.
            %       If they have 2, it is cyto, bg.
            if ch==1
                % Look for NADH first
                if exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_N_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_N_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_F_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_F_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_P_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_P_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_U_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_U_Probabilities.tif']);
                else
                    error('No NADH, FAD, or User prob maps present.');
                end
            elseif ch==2
                % Look for FAD first
                if exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_F_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_F_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_N_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_N_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_P_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_P_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_U_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_U_Probabilities.tif']);
                else
                    error('No NADH, FAD, or User prob maps present.');
                end
            elseif ch==3
                if exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_U_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_U_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_N_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_N_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_F_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_F_Probabilities.tif']);
                elseif exist([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_P_Probabilities.tif'],'file')==2
                    Prob = imread([Probpaths{pos} 'fr' num2str(frames(i),'%05i') '_P_Probabilities.tif']);
                else
                    error('No NADH, FAD, or User prob maps present.');
                end
            end
            
            % Build mask structure from Ilastik probs
            jointmask = boolean(zeros(size(Prob,1),size(Prob,2)));
            % Figure out which one is background, use it to find boundaries.
            % That will have the larges number of total pixels (largest area).
            % (Later only probs within this mask will be used to extract decays).
            %             Tots = squeeze(sum(sum(Prob,2),1)); % Old method of looking for largest area, but fails when sample is large, like CC's
            %             BG2 = Prob(:,:,Tots==max(Tots));
            BG = Prob(:,:,3)+Prob(:,:,4); % Exclude background and 'other' catmskory
            jointmask(BG<.3*255)=1; % Err on the side of excluding nucleus or BG
            jointmask = boolean(imfilter(jointmask,G,'same')); % Gauss blur
            jointmask(CRPind{ch})=0;
            % Load jointmasks into a cell
            jointmasks{frames(i),ch} = bwareafilt(jointmask,area_cuts);
            
            % Remove any formskround from BG mask.
            %             bgmasks{ch}(find(jointmasks{frames(i),ch})) = 0;
            bgmasks{ch}(find(jointmask)) = 0;
            
            %Get a rough mask of mito and cyto, joint, to autosclae the
            %images to something consistent (robust against bright crap).
            JntMsk = boolean(zeros(size(Prob,1),size(Prob,2)));
            Jnt = Prob(:,:,1)+Prob(:,:,2);
            JntMsk(Jnt>.7*255)=1;
            JntMsk = boolean(imfilter(JntMsk,G,'same')); % Gauss blur
            JntMsk(CRPind{ch})=0;
            JntMsks{i,ch} = bwareafilt(JntMsk,area_cuts);
            
            % Sort into labeled blobs
            % If merging blobs from a previous run, construct 'masks' cell from
            % those previous masks. Otherwise, get 'masks' from 'jointmasks{frames(i)}'
            % NOTE: this is currently only used for high-throughput dev
            % curves, which is just using NADH. Only code NADH case. Code
            % separate masks case if necessary, but it's complicated.
            % Separate masks case is mostly going to be needed for z-scans
            % of blasts, which will use SingleMasksBool. So just set
            % MergeBool=0 for these cases.
            if (exist([acqpath 'mask_merges.xls'])==2)&MergeBool
                clear Blobs2Merge masks
                if ~isempty(find(MergeFrs==frames(i)))
                    % Initialize masks to PrevMasks.
                    for blob = 1:length(PrevMasks{frames(i),ch}) masks{blob} = PrevMasks{frames(i),ch}(blob).L; end
                    
                    % Figure out which groups of fragments need to be merged
                    MergeInds = find(MergeFrs==frames(i)); %Indices for merges in this frame
                    DelInds=[];
                    for j = 1:length(MergeInds) % loop over embryos, may be multiple fragmented in one frame
                        MergedBlob = zeros(ydim,xdim1ch);
                        Blobs2Merge = MergeArr{MergeInds(j)};
                        for blob = Blobs2Merge % loop over fragments
                            MergedBlob = MergedBlob + PrevMasks{frames(i),ch}(blob).L;
                        end
                        
                        % Find blob with lowest number, make that the new merged blob
                        blob1 = min(Blobs2Merge); Blobs2Merge(Blobs2Merge==min(Blobs2Merge))=[];
                        masks{blob1} = MergedBlob;
                        
                        % Keep track of remaining blob indices to delete
                        DelInds = [DelInds Blobs2Merge];
                    end
                    masks(DelInds)=[];
                else
                    masks = bw2MasksCell(jointmasks{frames(i),ch});
                end
            else
                masks = bw2MasksCell(jointmasks{frames(i),ch});
            end
            masks(cellfun('isempty',masks))=[];
            
            % Optional auto-merge feature. If # of embs increased since
            % last frame, find which one fractured and group them. Use
            % area overlap and CoM distances
            
            if i>1&~JntMaskBool&~isempty(Masks0{frames(i-1)})
                clear ovlap CoMd blobmatch masks_merged
                
                if ~JntMaskBool
                    % Check that same number of masks were found in NADH and
                    % FAD images
                    if Chind(1)==1 & Chind(2)==2
                        if length(Masks0{frames(i-1),1})~=length(Masks0{frames(i-1),2})
                            error('Different number of masks found in NADH and FAD images. Check channels or use combined masks')
                        end
                    end
                end
                
                for blob = 1:length(Masks0{frames(i-1)},1) LastFrMasks{blob} = Masks0{frames(i-1),ch}(blob).L; end
                if size(masks,2)>size(LastFrMasks,2)
                    for blob = 1:length(masks)
                        for lastblob = 1:length(LastFrMasks)
                            ovlap(blob,lastblob) = length(intersect(find(masks{blob}),find(LastFrMasks{lastblob})));
                            [Cy1,Cx1] = imCoM(masks{blob}); [Cy2,Cx2] = imCoM(LastFrMasks{lastblob});
                            CoMd(blob,lastblob) = sqrt((Cy2-Cy1)^2+(Cx2-Cx1)^2);
                        end
                    end
                    
                    % First try to connect fragments to previous masks using maximum
                    % overlap, since that's the most robust. Ie, if they are overlapping,
                    % 99% of the time, they should be the same mask.
                    % If no overlap (e.g. a new blast fragment rmskistered), then tie the
                    % fragment to the clostest CoM
                    for blob = 1:size(ovlap,1)
                        if ~isempty(find(ovlap(blob,:)))
                            blobmatch(blob) = find(ovlap(blob,:)==max(ovlap(blob,:)));
                        else
                            blobmatch(blob) = find(CoMd(blob,:)==min(CoMd(blob,:)));
                        end
                    end
                    % Finally build cell of merged masks
                    for blob = 1:length(LastFrMasks)
                        mergeinds = find(blobmatch==blob);
                        MergedBlob = zeros(ydim,xdim1ch);
                        for k = 1:length(mergeinds)
                            MergedBlob = MergedBlob + masks{mergeinds(k)};
                        end
                        masks_merged{blob} = MergedBlob;
                    end
                    masks = masks_merged;
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
            
            % FINAL MASKS. Plug them into the struct
            if ~isempty(masks)
                for msk = 1:size(masks,2)
                    Masks0{frames(i),ch}(msk).L = masks{msk};
                    Masks0{frames(i),ch}(msk).Lper = bwperim(masks{msk});
                    [y,x] = find(masks{msk});
                    Masks0{frames(i),ch}(msk).CoM(1,:) = [mean(x) mean(y)];
                    Masks0{frames(i),ch}(msk).DdCoM(1,:) = [mean(x) mean(y)];
                    Coords{ch} = [Coords{ch}; mean(x) mean(y) frames(i)];
                    
                    % Looks like very bright pixels, especially in FAD channel,
                    % might be indicative of apoptotic activity. So store where
                    % the bright rmskions are for now. Maybe later can do an
                    % analysis. Nuclear rmskions might also be useful for
                    % something. Note: same array for all cells for a given
                    % frame. No other obvious way to store them.
                    Masks0{frames(i),ch}(msk).AboveThresh = AboveThresh;
                    Masks0{frames(i),ch}(msk).BelowThresh = BelowThresh;
                    %                 Masks0{frames(i),ch}(msk).NucInd = nucind;
                end
            else
                Masks0{frames(i),ch} = [];
            end
            1;
        end
        
        
        %         % Special case, if only one channel found a mask, don't use that
        %         % frame... not as relevant for human embryo stuff because NADH
        %         % and FAD images are pretty different, especially with Z-scans
        %         if NADHBool&FADBool
        %             if isempty(Masks0{frames(i),1})|isempty(Masks0{frames(i),2})
        %                 for c = Chind Masks0{frames(i),c}=[]; end
        %             end
        %         end
        1;
    end
    
    % Calculate background intensity values over time, using joint BG mask
    ImAllT = double(ones(iminf.Height,iminf.Width));
    % Erode to make sure we're far from the formskround
    for ch=1:3 BGmasks{pos,ch} = imerode(bgmasks{ch},strel('disk',10)); end
    
    for i = 1:length(frames)
        ImAllT = ImAllT + ims{frames(i)};
        for ch = Chind
            if ch==1
                if NADHBool&FADBool IntIm = ims{frames(i)}(:,1:ydim); else IntIm = ims{frames(i)}; end
            elseif ch==2
                if NADHBool&FADBool IntIm = ims{frames(i)}(:,ydim+1:end); else IntIm = ims{frames(i)}; end
            elseif ch==3
                IntIm = ims{frames(i)};
            end
            numscans = double(nameinds{(PosInd)&strcmp(nameinds(:,4),ChLabs{ch})'&([nameinds{:,6}]==frames(i)),9});
            IntVals = IntIm(find(BGmasks{pos,ch}));
            IntValsPerScan = IntVals/numscans;
            % BGvals matrix dims: [ValType, framenum, pos, ch];
            %    -> ValTypes: 1=intensity, 2=std, 3=NumOfPix, 4=timestamp
            BGvals(1,frames(i),pos,ch) = mean(IntValsPerScan); %
            BGvals(2,frames(i),pos,ch) = std(IntValsPerScan); %
            BGvals(3,frames(i),pos,ch) = length(find(BGmasks{pos,ch})); % BGnum_pixels
            BGvals(4,frames(i),pos,ch) = double(nameinds{(PosInd)&strcmp(nameinds(:,4),ChLabs{ch})'&([nameinds{:,6}]==frames(i)),8}); %
        end
    end
    
    % Display BG and save
    h = figure;
    imshow(ImAllT,[0 max(ImAllT(:))/2]); hold on
    [y1,x1] = find(bwperim(BGmasks{pos,1})); [y2,x2] = find(bwperim(BGmasks{pos,2})); [y3,x3] = find(bwperim(BGmasks{pos,3}));
    if NADHBool&FADBool plot(x1,y1,'b.'); plot(x2+ydim,y2,'g.');
    elseif NADHBool & ~FADBool plot(x1,y1,'b.');
    elseif ~NADHBool & FADBool plot(x2,y2,'g.');
    else plot(x3,y3,'r.');
    end
    saveas(h,[acqpath 'BG_masks.png']); close(h)
    
    
    % TRACKING: Calculate the trajectories using Dan Blaire's simple code
    param.mem = 3; param.good = 1; param.dim = 2; param.quiet = 0;
    
    % Do tracking for each channel, then so some sorting to make sure NADH
    % and FAD masks correspond correctly.
    for ch = Chind
        if ~JntMaskBool & length(frames)>1
            trajs{ch} = trackBlair(Coords{ch},50,param);
        end
    end
    
    imh = figure;
    set(gcf,'PaperPositionMode', 'auto')
    
    % Final sorting and displaying of masks, and ROI image saving
    for i = 1:length(frames)
        for ch = Chind
            % Bad frames?
            if ~isempty(Masks0{frames(i),ch})
                % Resort the structure mask elements according to particle IDs.
                for msk = 1:size(Masks0{frames(i),ch},2)
                    if ~JntMaskBool & length(frames)>1
                        tmptrj = trajs{ch};
                        trajID = tmptrj(tmptrj(:,1)==Masks0{frames(i),ch}(msk).DdCoM(1,1)&...
                            tmptrj(:,2)==Masks0{frames(i),ch}(msk).DdCoM(1,2)&tmptrj(:,3)==frames(i),4);
                        Masks{frames(i),ch}(trajID) = Masks0{frames(i),ch}(msk);
                    else
                        Masks = Masks0;
                        trajs = -1;
                    end
                end
                % Make a cell of joint masks as well in case you want to
                % analyze all embryos together.
                JntMasks{frames(i),ch}(1).L = jointmasks{frames(i),ch};
                JntMasks{frames(i),ch}(1).Lper = bwperim(jointmasks{frames(i),ch});
                JntMasks{frames(i),ch}(1).AboveThres = Masks0{frames(i),ch}(1).AboveThresh;
                JntMasks{frames(i),ch}(1).BelowThresh = Masks0{frames(i),ch}(1).BelowThresh;
            else
                for msk = 1:size(Masks0{frames(i),ch},2)
                    Masks{frames(i),ch}(msk).L = [];
                    Masks{frames(i),ch}(msk).Lper = [];
                end
                nameinds((PosInd)&strcmp(nameinds(:,4),ChLabs{ch})'&([nameinds{:,6}]==frames(i)),[7 8]) = num2cell(-1);
            end
        end
        
        % Additional sort in case NADH and FAD ended up with different mask sorting.
        if NADHBool & FADBool & ~JntMaskBool
            for msk = 1:size(Masks{frames(i),1},2) % masks in NADH channel
                [y,x] = find(Masks{frames(i),1}(msk).L); Ncom = [mean(x) mean(y)];
                for msk2 = 1:size(Masks{frames(i),2},2) % masks in FAD channel
                    [y2,x2] = find(Masks{frames(i),2}(msk2).L); Fcom = [mean(x2) mean(y2)];
                    % Find nearest mask to current NADH mask
                    dists(msk2) = sqrt((Fcom(1)-Ncom(1))^2+(Fcom(2)-Ncom(2))^2);
                end
                minind = find(dists==min(dists));
                % If indices don't correspond, re-order FAD masks so they do
                if minind~=msk
                    tmp = Masks{frames(i),2}(msk).L; tmpper = Masks{frames(i),2}(msk).Lper;
                    Masks{frames(i),2}(msk).L = Masks{frames(i),2}(msk2).L;
                    Masks{frames(i),2}(msk).Lper = Masks{frames(i),2}(msk2).Lper;
                    Masks{frames(i),2}(msk2).L = tmp;
                    Masks{frames(i),2}(msk2).Lper = tmpper;
                end
            end
        end
        
    end
    
    % Finish plotting ROI images and save
    for i = 1:length(frames)
        xdim = size(ims{frames(i)},2); ydim = size(ims{frames(i)},1);
        % Get foreground intensities to do custom display sclaing
        clear Ints JntDbl
        if NADHBool&FADBool % Dual image
            JntDbl = [JntMsks{i,1} JntMsks{i,1}];
            Ints = ims{frames(i)}(find(JntDbl));
        else % Individual channels (e.g. only NADH, or FAD plus ruthenium)
            for ch = Chind
                Ints(:,ch) = ims{frames(i)}(find(JntMsks{i,ch}));
            end
            Ints = reshape(Ints,[],1);
        end
        lb = max([0 mean(Ints)-2*std(Ints)]); ub = min([200 (mean(Ints)+4*std(Ints))]);
        imshow(ims{frames(i)},[lb ub],'Border','tight'); hold on;
        if ~ProcDisp
            saveas(imh,[ROIpaths{pos} 'ROI_im' num2str(frames(i),'%03i') '.png'],'png');
        end
        
        % Dual channel (common) case
        if NADHBool&FADBool
            for ch = Chind
                % If no masks were found in this channel, show in ROI image
                if isempty(Masks{frames(i),ch})
                    text(xdim1ch*(ch-1)+xdim1ch/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',18)
                else
                    for msk = 1:size(Masks{frames(i),1},2)
                        [y,x] = find(Masks{frames(i),1}(msk).Lper);
                        Ints = ims{frames(i)}(find([JntMsks{i,1} JntMsks{i,2}]));
                        plot(x,y,'b.','markersize',3)
                        if ~JntMaskBool text(mean(x),mean(y),num2str(msk),'HorizontalAlignment','center','fontsize',14,'color','r'); end
                    end
                    for msk = 1:size(Masks{frames(i),2},2)
                        [y,x] = find(Masks{frames(i),2}(msk).Lper);
                        plot(x+xdim1ch,y,'g.','markersize',3)
                        if ~JntMaskBool text(mean(x)+xdim1ch,mean(y),num2str(msk),'HorizontalAlignment','center','fontsize',14,'color','r'); end
                    end
                end
            end
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            set(gcf,'paperPositionMode','auto')
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            set(gcf,'inverthardcopy','off')
            saveas(imh,[ROIpaths{pos} 'ROI_im' num2str(frames(i) ,'%03i') '.png'],'png');
            hold off;
        else
            for ch = Chind
                if isempty(Masks{frames(i),ch})
                    text(xdim1ch*(ch-1)+xdim1ch/2,ydim/2,'No mask found','color','red','HorizontalAlignment','center','fontsize',18)
                else
                    for msk = 1:size(Masks{frames(i),ch},2)
                        [y,x] = find(Masks{frames(i),ch}(msk).Lper);
                        plot(x,y,'r.','markersize',3)
                        if ~JntMaskBool text(mean(x),mean(y),num2str(msk),'HorizontalAlignment','center','fontsize',12,'color','r'); end
                    end
                end
            end
            set(gcf,'paperPositionMode','auto')
            set(gcf,'position',[200 200 xdim*.6 ydim*.6]);
            set(gcf,'inverthardcopy','off')
            saveas(imh,[ROIpaths{pos} 'ROI_im' num2str(frames(i),'%03i') '_' ChLabs{ch} '.png'],'png');
            hold off;
        end
    end
    
    
    %     else % whole frame ROIs... finish and test at some point
    %         for i = 1:length(frames)
    %             for ch = Chind
    %                 Masks{frames(i),ch}(1).L = ones(ydim,ydim);
    %             end
    %         end
    %     end
    
    
    if ~JntMaskBool save([acqpath 'Masks_' Dpos(pos).name '.mat'],'frames','Masks','Coords','trajs'); end
    Masks = JntMasks;
    save([acqpath 'Joint' MaskLab '_' Dpos(pos).name '.mat'],'frames','Masks','Coords','trajs');
end

save([acqpath 'multiD_indices.mat'],'nameinds')
save([acqpath 'BG_vals.mat'],'BGvals')
% save([acqpath 'BG_masks_ints.mat'],'BGmasks','BGirrs','BGstds','BGnum_pixels','BGtimestps')

% To plot BG intensities, use 'BG_ints_plot(acqpath,pos)'
if NADHBool&FADBool BG_ints_plot(acqpath); end % Only coded for dual case, not worth adapting to single.


