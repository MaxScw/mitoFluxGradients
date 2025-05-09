function FLIM_decay_from_probmap(acqpath,poss,MaskLab,MinFrames,ProbThreshMeth,MasksToDo)
% Extracts FLIM decay curves from FLIMages and probability maps obtained
% with Ilastik. 
% INPUTs:
% -acqpath: path to acquisition
% -poss: path to acquisition
% -MaskLab: path to acquisition
% -MinFrames: path to acquisition
% -ProbThreshMeth - 1=prob map, 2=prob^2, any number between 0 and 1
%   specifies a custom prob threshold. DEFAULT value is 0.7 hard thresh.
% -MasksToDo - cell of arrays of mask numbers to get decays for. Each cell
%   element has an array of mask numbers for the corresponding position. For
%   a given position, if you want to do all the masks, enter '-1'.
%   NOTE: array must be in the same order as the '..._Masks.mat' files in the
%   acquisition acqpath. For example, if there are Pos0, Pos1, and Pos2, but
%   Pos0 didn't have any data, you would order a 2-element cell corresponding
%   to {[Pos1MasksToDo],[Pos2MasksToDo]}

% % TEST in script mode:
% clear all;
% acqpath = 'C:\Users\Tim\Documents\Academic - Research\Publication Materials\2019 Denny Intro to MetIm\Figures\Mitotracker\Emb1_a3_5nM_dual_20s_42MitoHWP_26FADHWP_z0';
% MaskLab = 'JointMasks';
% MasksToDo = {[-1],[2 4],[],[],[],[],[]};

if acqpath(end)~='\' acqpath = [acqpath '\']; end;
slashes = strfind(acqpath,'\');
Run = acqpath(slashes(end-1)+1:end-1);
% New default for Boston IVF, typically imaging 1 
if ~exist('MaskLab')|MaskLab==-1 MaskLab = 'JointMasks'; end 
% if ~exist('MaskLab')|MaskLab==-1 MaskLab = 'mask'; end
if ~exist('MinFrames')|MinFrames==-1 MinFrames = 1; end
if ~exist('ProbThreshMeth')|ProbThreshMeth==-1 ProbThreshMeth = .7; end
IRfact=1;

% Trying to be too clever here, I think
% % Check if only 'JointMasks' were produced, then assume MaskLab = 'JointMasks'
% Dmsk = dir([acqpath 'Masks*.mat']); Djnt = dir([acqpath 'JointMasks_*.mat']);
% if isempty(Dmsk)&~isempty(Djnt) MaskLab = 'JointMasks'; end

% Load nameinds to be sure to connect correct mask to each sdt file
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end

% Determine which channels are present
uchans = unique(nameinds(:,4));
if ~isempty(find(strcmp(uchans,'NADH'))) Chind(1) = 1; end
if ~isempty(find(strcmp(uchans,'FAD'))) Chind(2) = 2; end
if ~isempty(find(strcmp(uchans,'UserChan'))) Chind(3) = 3; end
Chind(Chind==0)=[];
ChLabs = {'NADH','FAD','UserChan'};


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

sdtpath = [acqpath 'sorted_sdts\'];
Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    % Check that folder is a 'Pos', and see if there is a corresponding 'Mask' file
    if ~strcmp(Dpos(i).name(1:3),'Pos')|isempty(dir([acqpath MaskLab '_' Dpos(i).name  '.mat']))
        remove = [remove i];
    end
end
Dpos(remove)=[];

for pos = 1:size(Dpos,1)
    
    clear decay_structs
    PosNum = Dpos(pos).name(4:end);
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),PosNum)
            continue;
        end
    end
    
    % Load Masks
    load([acqpath MaskLab '_' Dpos(pos).name  '.mat']);
    try
        WowCropMask = tiffread2([acqpath 'sorted_sdts\Crop_' Dpos(pos).name '_' Run '.tif']);
        WowCropMask = boolean(WowCropMask.data);
    catch
        disp('No ImageJ Crop file present')
    end
    
    %% Load all data
    Dsdt = dir([acqpath 'sorted_sdts\' Dpos(pos).name '\*.sdt']);
    NumOfStds = length(Dsdt);
    filenames = [];
    clear dashes t chans z StkFr
    for i = 1:NumOfStds
        filenames{i}=Dsdt(i).name;
        dashes = strfind(filenames{i},'_');
        t(i) = str2num(filenames{i}(dashes(1)+1:dashes(2)-1));
        chans{i} = filenames{i}(dashes(2)+1:dashes(3)-1);
        chnums(i) = find(strcmp(ChLabs,chans{i}));
        % Will be the same indexing for all msks
        z(i) = str2num(filenames{i}(dashes(3)+1:dashes(3)+4));
        
        % Use these and nameinds to get the stack slice number, which
        % corresponds with the 'Masks' structure elements, used below
        subset = find((PosInd)&([nameinds{:,2}]==t(i))&([nameinds{:,5}]==z(i))&([nameinds{:,7}]>-1));
        
        % For NADH and FAD, there should be only 2 frames. If 1 or less, a 
        % frame got dropped, so just set decay_structs{i,msk} to [];
        if isempty(strfind([nameinds{subset,4}],chans{i}))
            StkFr(i) = -1; % flag
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
    nummsks = 0;
    for i = 1:size(Masks,1) nummsks(i) = size(Masks{i,chnums(i)},2); end
    nummsks = max(nummsks);
    decay_structs = cell(NumOfStds,nummsks);
    
    %Two photon image block
    block=1; %1:2pf, 2:SHG 
    % NOTE: If you are ever doing something like 2-color FLIM with a
    % dichroic and need to FLIM analyze the data in ch2, temporarily change
    % block to 2, but make sure to change it back!
    pathpar = [acqpath 'sorted_sdts\' Dpos(pos).name '\'];
    filenames = filenames;
    for i = 1:NumOfStds
        % use the StkFr, -1 flag here to skip
        % frames that dropped a frame or had something go wrong.
        if StkFr(i)<0
            continue;
        end
        disp([Run ', Pos' PosNum ', Frame ' num2str(StkFr(i)) ', ' ChLabs{chnums(i)}]);
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
        %         acqpath = acqpath;
        dt = dt;
        %image plot handle
        image_handle = -99;
        %1 if pixel selected, 0 if not
        selected_pixel = zeros(size(img));
        %Handles for plot showing selected pixels
        selected_pixel_handle = [];
        %FLIM Data
        flim = uint8(flim);
        % Fluorescence decay data extracted from selected pixels
        % We will be segmenting mitochondria and cytoplasm now, so need
        % separate decays. Still keep the ability to look at cyto and mito
        % together using 'decayjoint'.
        % NOTE: for low res where you only have foreground and background,
        % foreground will be the 'cyto' channel
        decaymito = zeros(size(flim,1),1);
        decaycyto = zeros(size(flim,1),1);
        decayjoint = zeros(size(flim,1),1);
        
        
        %% Get decay for image from pixels contained in the mask.
        
        % Do for each msk. For single generic mask for whole pictures (eg
        % cumulus cells), store mask in first element of 'Masks' cell.
        if exist('MasksToDo')&(MasksToDo{pos}~=-1)
            MaskInd = MasksToDo{pos};
        else
            MaskInd = 1:size(Masks{StkFr(i),chnums(i)},2);
        end
        
        % Load the probability maps corresponding to this sdt.
        ProbPath = [sdtpath 'ProbMaps_Pos' PosNum '_' Run '\'];
        if exist(ProbPath) % Make sure prob maps were calculated
            if strcmp(chans{i},'FAD')
                % Try to load FAD prob map, but use NADH one if that's all
                % that's there
                if exist([ProbPath 'fr' num2str(StkFr(i),'%05i') '_F_Probabilities.tif'])==2
                    Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_F_Probabilities.tif']);
                elseif exist([ProbPath 'fr' num2str(StkFr(i),'%05i') '_N_Probabilities.tif'])==2
                    Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_N_Probabilities.tif']);
                elseif exist([ProbPath 'fr' num2str(StkFr(i),'%05i') '_P_Probabilities.tif'])==2
                    Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_P_Probabilities.tif']);
                else
                    error('Channel problem')
                end
            elseif strcmp(chans{i},'NADH')
                if exist([ProbPath 'fr' num2str(StkFr(i),'%05i') '_N_Probabilities.tif'])==2
                    Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_N_Probabilities.tif']);
                elseif exist([ProbPath 'fr' num2str(StkFr(i),'%05i') '_F_Probabilities.tif'])==2
                    Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_F_Probabilities.tif']);
                elseif exist([ProbPath 'fr' num2str(StkFr(i),'%05i') '_P_Probabilities.tif'])==2
                    Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_P_Probabilities.tif']);
                else
                    error('Channel problem')
                end
            elseif strcmp(chans{i},'UserChan')
                Prob = imread([ProbPath 'fr' num2str(StkFr(i),'%05i') '_U_Probabilities.tif']);
            else
                error('Problem with channel masks')
            end
        else % If only binary masks were created, set Prob to be equal to the masks
            Prob = zeros([size(Masks{StkFr(i),chnums(i)}(1).L) 3]);
            for msk = 1:size(Masks{StkFr(i),chnums(i)},2) Prob(:,:,1) = Prob(:,:,1) + Masks{StkFr(i),chnums(i)}(msk).FL; end
            Prob = Prob.*255;
        end
        
        for msk = MaskInd
            %time axis
            decay = zeros(size(flim,1),1);
            time = (1:length(decay))'*dt;
            
            % If this msk wasn't found for this frame, continue
            if msk>size(Masks{StkFr(i),chnums(i)},2)
                continue;
            end
            fn = fieldnames(Masks{StkFr(i),chnums(i)}(msk));
            if isempty(eval(['Masks{StkFr(i),chnums(i)}(msk).' fn{1}]))
                continue;
            end
            
            % Load mask into selected_pixel, used to extract FLIM decay
            selected_pixel=Masks{StkFr(i),chnums(i)}(msk).L;
            
            % DECAY EXTRACTION
            % If pixels found in mask
            if (~isempty(selected_pixel))&(StkFr(i)~=-1)
                % WEIGHTING ARRIVAL TIMES: For each channel (mito, cyto), use
                % the probs within the mask for this msk, and weight the
                % arrival times according to the probs of that pixel.
                % Make 3D matrices of repeated prob maps and masks, then
                % simply multiply them by the FLIMage to get a prob-weighted
                % decay from only within that mask.
                selected_pixel = double(selected_pixel);
                ProbMito = double(squeeze(Prob(:,:,1)));
                ProbCyto = double(squeeze(Prob(:,:,2)));
                
                % Check for ProbBestExps file, containing exponents for prob map adjustment
                if exist([acqpath 'ProbBestExps_Pos' PosNum '.mat'])
                    load([acqpath 'ProbBestExps_Pos' PosNum '.mat'])
                    MitoExp = BestExps(1,StkFr(i));
                    CytoExp = BestExps(2,StkFr(i));
                    ProbMito = ProbMito.^MitoExp; ProbMito = ProbMito./max(ProbMito(:)).*255;
                    ProbCyto = ProbCyto.^CytoExp; ProbCyto = ProbCyto./max(ProbCyto(:)).*255;
                end
                
                if ProbThreshMeth==1
                        % Prob method. No adjustment necessary. Lifetimes from each 
                        % pixel are weighted according to their probability of 
                        % being 'mitochondria', 'cytoplasm', or 'background'
                        % Note: Tends toward higher cross-contamination
                        decay_structs{i,msk}.SegMeth = 'ProbWeight';
                elseif ProbThreshMeth==2 % Tried. Not using.
                        % Prob^2 method - sharpen boundaries between segments
                        ProbMito = ProbMito.^2./255;
                        ProbCyto = ProbCyto.^2./255;
                        decay_structs{i,msk}.SegMeth = 'Prob^2Weight';
                elseif ProbThreshMeth>0&ProbThreshMeth<1
                        % Custom hard thresh method. 
                        ProbMito(ProbMito>=ProbThreshMeth*255)=255; ProbMito(ProbMito<ProbThreshMeth*255)=0;
                        ProbCyto(ProbCyto>=ProbThreshMeth*255)=255; ProbCyto(ProbCyto<ProbThreshMeth*255)=0;
                        decay_structs{i,msk}.SegMeth = '0p7 prob thresh';
                else
                    error('Problem with your prob thresh method.')
                end
                for at=1:size(flim,1)
                    decaymito(at) = sum(sum(selected_pixel.*ProbMito.*double(squeeze(flim(at,:,:)))))/255;
                    decaycyto(at) = sum(sum(selected_pixel.*ProbCyto.*double(squeeze(flim(at,:,:)))))/255;
                    % divide by 255 to convert to probs, hence covert
                    % decays to units of estimated # of photons
                end
                decayjoint = decaymito + decaycyto;
                decay = [decayjoint decaymito decaycyto]; % Keep them all in this one variable as separate columns
                % NOTE: here we are switching the order so that 'joint'
                % comes first. For single segment data, this will be the
                % only segment.
                
                % Calculate irradiance and irradiance std for averages
                % later. Get intensity at each pixel, I(x,y). Then get
                % the average for each class by weighting that
                % intensity by the prob.
                ind = find(selected_pixel);
                % Scale by IllProf
                IntIm = squeeze(sum(flim,1))./IllProfCal.*mean(IllProfCal(:));
                IntVals = IntIm(ind);
                IntValsPerScan = IntVals/numscans;
                ProbValsMito = double(ProbMito(ind))./255;
                ProbValsCyto = double(ProbCyto(ind))./255;
                photons = sum(IntVals);
                irr(1) = (sum(IntValsPerScan.*ProbValsMito)+sum(IntValsPerScan.*ProbValsCyto))/(sum(ProbValsMito)+sum(ProbValsCyto));
                irr(2) = sum(IntValsPerScan.*ProbValsMito)/sum(ProbValsMito);
                irr(3) = sum(IntValsPerScan.*ProbValsCyto)/sum(ProbValsCyto);
                irr_std(1) = sqrt(var(IntValsPerScan,ProbValsMito+ProbValsCyto));
                irr_std(2) = sqrt(var(IntValsPerScan,ProbValsMito));
                irr_std(3) = sqrt(var(IntValsPerScan,ProbValsCyto));
                
                %update total photon counts
                totphot = num2str(sum(decay(:,3)));
                
                %% Save decays in the same format that is output by 'FLIMDataAnalysisGUI_ver5_3'
                decay_structs{i,msk}.decay = decay; % Keep 'decay' in there for backwards compatibility with other functions that look for it.
                decay_structs{i,msk}.name = filename;
                decay_structs{i,msk}.filename = filename;
                decay_structs{i,msk}.image = image;
                decay_structs{i,msk}.selected_pixel = selected_pixel;
                decay_structs{i,msk}.time = time;
                decay_structs{i,msk}.num_pixel = sum(sum(selected_pixel));
                decay_structs{i,msk}.numscans = numscans;
                % CalInt previously used to keep track of laser output
                % power and adjust the measured intensities accordingly.
                % Stopped using it, so just set to -1 for now, possible
                % revive in the future to make intensities as accurate as
                % possible.
                decay_structs{i,msk}.CalInt = -1; 
                decay_structs{i,msk}.photons = photons;
                decay_structs{i,msk}.irr = irr;
                decay_structs{i,msk}.irr_std = irr_std;
                decay_structs{i,msk}.IRfact = IRfact;
                decay_structs{i,msk}.timestp = timestp(i);
                decay_structs{i,msk}.segments = 'joint, mito, cyto';
                decay_structs{i,msk}.fit_region  = [round(resol*(1-meas.tac_lh/100)+1) round(resol*(1-meas.tac_ll/100)-1)];
                decay_structs{i,msk}.noise_region = [round(resol*(1-meas.tac_lh/100)+1) round(resol*(1-meas.tac_lh/100)+1)+5];
            end
        end
        fclose('all');
    end
    
    
    % Final save for this position
    for msk = MaskInd
        decay_struct = cell(size(nameinds,1),1);
        %         decay_struct = decay_structs(:,msk);
        % Index decays identically to column 1 of nameinds, with blanks for
        % missing decays
        for i = 1:NumOfStds
            decaysind = PosInd&[nameinds{:,2}]==t(i)&[nameinds{:,5}]==z(i)&strcmp(nameinds(:,4),chans{i})';
            decay_struct(decaysind) = decay_structs(i,msk);
        end
        
        decay_struct
        
        
        % Filter mask numbers that were short lived, probably false mask
        % segmentation
        NonemptyFrms = cell2mat(nameinds(find(~cellfun('isempty',decay_struct)),6));
        if length(unique(NonemptyFrms))>=MinFrames
            % Save decays
            if strcmp(MaskLab,'Masks')
                save([acqpath 'decays_' Dpos(pos).name '_mask' num2str(msk) '.mat'],'decay_struct');
            else
                save([acqpath 'decays_' Dpos(pos).name '_' MaskLab '.mat'],'decay_struct');
            end
        end
    end
end
