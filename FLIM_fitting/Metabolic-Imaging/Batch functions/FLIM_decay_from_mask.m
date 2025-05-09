function FLIM_decay_from_mask(path,poss,MaskLab,MinFrames,MasksToDo)
% Performs the same functions as 'FLIMDataAnalysisGUI_ver5_3', but without
% the GUI, and with preset params so that it can be operated in batch mode,
% so I can go have a beer while it chimes out the results.
% This version takes a mask of which pixels are to be binned. Mask is
% produced from 'EllMask_Iner', which does a bpass, threshold, to find
% roughly circular blobs. Thus, no threshold is necessary.
% The input is the mask elements in linear index form (like, raster in a 2d
% image, like the 'find' matlab function)
% Meant to run on semi-isolated eggs/embryos. If blobs are overlapping, it
% won't really work.

% clear all;
% path = 'C:\Users\Tim\Documents\Academic - Research\Data\2017-07-08 Batch 8\s1_a4_postDCFDA\';
% MaskLab = 'SingleMasks';
% MasksToDo = {[-1],[2 4],[],[],[],[],[]};

% MasksToDo - cell of arrays of mask numbers to get decays for. Each cell
% element has an array of mask numbers for the corresponding position. For
% a given position, if you want to do all the masks, enter '-1'.
% NOTE: array must be in the same order as the '..._Masks.mat' files in the
% acquisition path. For example, if there are Pos0, Pos1, and Pos2, but
% Pos0 didn't have any data, you would order a 2-element cell corresponding
% to {[Pos1MasksToDo],[Pos2MasksToDo]}

% Versions
% 2018-08-03: Moved to probability framework, but had to update this to
% reanalyze Emre data with corrected TAC and dt_irf
% 2017-07-14: Add a condition to handle arbitrary user channel (not NADH or
%  FAD)
% 2017-02-22: Add a conditional to look for channels with ruthenium dye. If
% they are there, calculate fluorescence intensity in the ruthenium
% channel, as well as the FAD channel in the extracellular areas. 
% 2015-12-29: Got macro to start logging IRpow from dumbass Spectra MaiTai
% GUI. Saves in a text file, then can use that additional info to scale
% irradiances in this routine.
% 2015-09-23: Changed 'Mask' structure in 'Make_Masks_Eggs_DBTrack'. Had to
% change this to look for proper Mask structure.
% 2015-12-04: Adapt to include illumination profile calibration image.
% Scale intensity images to get scaled irradiance
% 2014-11-07 Rearrange so that parfors don't pass 500mb structures to each
% of the workers. Parfors are freezing.... Come back to this some time.
% 2014-10-10 Reindex decay_struct to match first column of nameinds


if path(end)~='\' path = [path '\']; end;
slashes = strfind(path,'\');
Run = path(slashes(end-1)+1:end-1);
if ~exist('MaskLab')|MaskLab==-1 MaskLab = 'Masks'; end
if ~exist('MinFrames')|MinFrames==-1 MinFrames = 1; end

% Load nameinds to be sure to connect correct mask to each sdt file
load([path 'multiD_indices.mat']);

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

sdtpath = [path 'sorted_sdts\'];
Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    if ~strcmp(Dpos(i).name(1:3),'Pos')
        remove = [remove i];
    end
    % See if there is a corresponding 'Mask' file
    if isempty(dir([path Run '_' Dpos(i).name '_' MaskLab '.mat']))
        remove = [remove i];
    end
end
Dpos(remove)=[];

for posnum = 1:size(Dpos,1)
    
    clear decay_structs
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
    
    % Load Masks
    load([path Run '_' Dpos(posnum).name '_' MaskLab '.mat']);
    try 
        WowCropMask = tiffread2([path 'sorted_sdts\Crop_' Dpos(posnum).name '_' Run '.tif']);
        WowCropMask = boolean(WowCropMask.data);
    catch
        disp('No ImageJ Crop file present')
    end
    
    %% 'Open Image' button. Load all data
    Dsdt = dir([path 'sorted_sdts\' Dpos(posnum).name '\*.sdt']);
    NumOfStds = length(Dsdt);
    filenames = [];
    clear dashes t chans z StkFr
    for i = 1:NumOfStds
        filenames{i}=Dsdt(i).name;
        dashes = strfind(filenames{i},'_');
        t(i) = str2num(filenames{i}(dashes(1)+1:dashes(2)-1));
        chans{i} = filenames{i}(dashes(2)+1:dashes(3)-1);
        % Will be the same indexing for all eggs
        z(i) = str2num(filenames{i}(dashes(3)+1:dashes(3)+4));
        
        % Use these and nameinds to get the stack slice number, which
        % corresponds with the 'Masks' structure elements, used below
        subset = find((PosInd)&([nameinds{:,2}]==t(i))&([nameinds{:,5}]==z(i))&([nameinds{:,7}]>-1));
        % Should be only 2 frames. If 1 or less, a frame got dropped, so
        % just set decay_structs{i,egg} to [];
        if isempty(strfind([nameinds{subset,4}],chans{i}))
            StkFr(i) = -1; % flag
        elseif strcmp(chans{i},'Ruth')
            StkFr(i) = nameinds{i,6};
            timestp(i) = nameinds{i,8};
        else
            if strfind(nameinds{subset(1),4},chans{i})
                StkFr(i) = nameinds{subset(1),6};
                timestp(i) = nameinds{subset(1),8};
            elseif strfind(nameinds{subset(2),4},chans{i})
                StkFr(i) = nameinds{subset(2),6};
                timestp(i) = nameinds{subset(2),8};
            else
                error('Something wrong with the name indexes')
            end
        end
    end
    uniqchans = unique(chans);
    
    %Image structure that contains all the information about the newly
    %opened images
    numeggs = 0;
    for i = 1:size(Masks,2) numeggs(i) = size(Masks{i},2); end
    numeggs = max(numeggs);
    decay_structs = cell(NumOfStds,numeggs);
    
    %Two photon image block
    block=1; %1:2pf, 2:SHG
    pathpar = [path 'sorted_sdts\' Dpos(posnum).name '\'];
    filenames = filenames;
    for i = 1:NumOfStds
        % use the StkFr, -1 flag here to skip
        % frames that dropped a frame or had something go wrong.
        if StkFr(i)<0
            continue;
        end
              %if StkFr(i)==2
               %     disp('')
              %end
        disp([Run ', Pos' uManPos ', Frame ' num2str(StkFr(i))]);
        %load sdt file
        sdt = bh_readsetup([pathpar filenames{i}]);
        AcqRng = GetPhotonCollectTRange([pathpar filenames{i}]);
        
        % Get indexes
        dashes = strfind(filenames{i},'_');
        
        ch = bh_getdatablock_v095(sdt,block);
        img = uint8(squeeze(sum(ch,1)));
        flim = ch;
        
        % How many scans were integrated to get total intensity for this frame
        % This is a stupid BH thing. Num of scans varies from frame to frame
        meas = bh_getmeasdesc(sdt,1);
        numscans = double(meas.hist_fida_points);
        
        %time/channel = range/(gain*ADCresolution)
        range = sdt.SP_TAC_R*10^9;
        gain = double(sdt.SP_TAC_G);
        resol = double(sdt.SP_ADC_RE);
        dt = range/(gain*resol);
        
        %%
        %Update image structure
        image = img;
        filename = filenames{i};
        %         path = path;
        dt = dt;
        %image plot handle
        image_handle = -99;
        %1 if pixel selected, 0 if not
        selected_pixel = zeros(size(img));
        %Handles for plot showing selected pixels
        selected_pixel_handle = [];
        %FLIM Data
        flim = flim;
        %Fluorescence decay data extracted from selected pixels
        decay = zeros(size(flim,1),1);
        
        
        %% Get decay for image from pixels contained in the mask.
        
        % Do for each egg. For single generic mask for whole pictures (eg
        % cumulus cells), store mask in first element of 'Masks' cell
        if exist('MasksToDo')&(MasksToDo{posnum}~=-1)
            MaskInd = MasksToDo{posnum};
        else
            MaskInd = 1:size(Masks{StkFr(i)},2);
        end
        for egg = MaskInd
            %time axis
            decay = zeros(size(flim,1),1);
            time = (1:length(decay))'*dt;
            
            % If this egg wasn't found for this frame, continue
            if egg>size(Masks{StkFr(i)},2)
                continue;
            end
            fn = fieldnames(Masks{StkFr(i)}(egg));
            if isempty(eval(['Masks{StkFr(i)}(egg).' fn{1}]))
                continue;
            end
            
            % Masks may have only one set of masks for both channels (FAD
            % and NADH). C elegans don't move, so this works well. Mouse
            % eggs move a bit, so it's better to find masks in respective
            % channels. Check for what time of mask is here and load
            % appropriate mask
            if isfield(Masks{StkFr(i)}(egg),'L')
                selected_pixel=Masks{StkFr(i)}(egg).L;
            else
                if strcmp(chans{i},'FAD')
                    selected_pixel=Masks{StkFr(i)}(egg).FL;
                elseif strcmp(chans{i},'NADH')
                    selected_pixel=Masks{StkFr(i)}(egg).NL;
                elseif strcmp(chans{i},'UserChan')
                    selected_pixel=Masks{StkFr(i)}(egg).NL;
                elseif strcmp(chans{i},'Ruth')
                    selected_pixel=imerode(WowCropMask&~Masks{StkFr(i)}(egg).NL,strel('disk',10));
                else
                    error('Problem with channel masks')
                end
                CalInt = -1;
            end
            
            % DECAY EXTRACTION
            % If pixels found in mask
            if (~isempty(selected_pixel))&(StkFr(i)~=-1)
                reshapedselectedpix = repmat(reshape(uint16(selected_pixel),[1,size(selected_pixel)]),[length(decay),1,1]);
                decay = decay + sum(sum(reshapedselectedpix.*flim,2),3);
                %                 close all;imshowt(selected_pixel); figure;imshowt(image); pause(.3)
                
                % Calculate irradiance and irradiance std for averages later
                % Note:
                % I previously had a 'ProfBool' here and 'irrSc' to be
                % 'scaled irradiance'. But now we always use illumination
                % profiles and the scaled irradiance. Unscaled irradiance
                % is not really meaningful or accurate, so for simplicity 
                ind = find(selected_pixel);
                IntIm = squeeze(sum(flim,1));
                % Scale by IllProf
                IntIm = IntIm./IllProfCal.*mean(IllProfCal(:));
                IntVals = IntIm(ind);
                photons = sum(IntVals);
                irr = mean(IntVals)/numscans;
                irr_std = std(IntVals)/numscans;
                
                %update total photon counts
                totphot = num2str(sum(decay));
                
                %% Save decays in the same format that is output by 'FLIMDataAnalysisGUI_ver5_3'
                decay_structs{i,egg}.decay = decay;
                decay_structs{i,egg}.name = filename;
                decay_structs{i,egg}.filename = filename;
                decay_structs{i,egg}.image = image;
                decay_structs{i,egg}.selected_pixel = selected_pixel;
                decay_structs{i,egg}.time = time;
                decay_structs{i,egg}.num_pixel = sum(sum(selected_pixel));
                decay_structs{i,egg}.numscans = numscans;
                decay_structs{i,egg}.CalInt = CalInt;
                decay_structs{i,egg}.photons = photons;
                decay_structs{i,egg}.irr = irr;
                decay_structs{i,egg}.irr_std = irr_std;
                decay_structs{i,egg}.timestp = timestp(i);
            end
        end
        fclose('all');
    end
    
    
    % Final save for this position
    for egg = MaskInd
        decay_struct = cell(size(nameinds,1),1);
        %         decay_struct = decay_structs(:,egg);
        % Index decays identically to column 1 of nameinds, with blanks for
        % missing decays
        for i = 1:NumOfStds
            decaysind = PosInd&[nameinds{:,2}]==t(i)&[nameinds{:,5}]==z(i)&strcmp(nameinds(:,4),chans{i})';
            decay_struct(decaysind) = decay_structs(i,egg);
        end
        
        % Filter mask numbers that were short lived, probably false mask
        % segmentation
        NonemptyFrms = cell2mat(nameinds(find(~cellfun('isempty',decay_struct)),6));
        if length(unique(NonemptyFrms))>=MinFrames
            % Save decays
            if strcmp(MaskLab,'Masks')
                save([path 'decays_' Dpos(posnum).name '_mask' num2str(egg) '.mat'],'decay_struct');
            else
                save([path 'decays_' Dpos(posnum).name '_' MaskLab '.mat'],'decay_struct');
            end
            
        end
    end
end
