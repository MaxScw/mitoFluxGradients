function FLIM_decay_from_mask_AnyFLIMageSeries(path,poss,MaskLab,MinFrames,MasksToDo)
% Generalization of 'FLIM_decay_from_mask.m' to analyze any arbitrary
% series of FLIMages in a folder ('path').
% IMPORTANT: You must run 'IntensitymatrixTiffs.m' and
% 'Make_Mask_Generic_WoW_AnyFLIMageSeries.m' on that path first to convert
% to tiffs and create masks.
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
% path = 'C:\Dropbox\data\2016-11-17ROS_Test_marta\better_ones\';
% MasksToDo = {[-1],[2 4],[],[],[],[],[]};

% MasksToDo - cell of arrays of mask numbers to get decays for. Each cell
% element has an array of mask numbers for the corresponding position. For
% a given position, if you want to do all the masks, enter '-1'.
% NOTE: array must be in the same order as the '..._Masks.mat' files in the
% acquisition path. For example, if there are Pos0, Pos1, and Pos2, but
% Pos0 didn't have any data, you would order a 2-element cell corresponding
% to {[Pos1MasksToDo],[Pos2MasksToDo]}

% Versions

if path(end)~='\' path = [path '\']; end;
slashes = strfind(path,'\');
Run = path(slashes(end-1)+1:end-1);

% Check for IllProfCal.mat file
Dprof = dir([UpOneDir(path) '\DailyFiles\IllProfCal.mat']);
ProfBool = 0;
if ~isempty(Dprof)
    load([UpOneDir(path) '\DailyFiles\IllProfCal.mat']);
    IllProfCal = double(IllProfCal);
    ProfBool = 1;
end

load([path Run '_SingleMasks.mat']);

Dsdt = dir([path '\*.sdt']);
NumOfStds = length(Dsdt);
filenames = [];
clear dashes t chans z StkFr

%Image structure that contains all the information about the newly
%opened images
numeggs = 1;
decay_structs = cell(NumOfStds,numeggs);

%Two photon image block
block=1; %1:2pf, 2:SHG
for i = 1:NumOfStds
    % Skip frames for which no mask could be found
    if isempty(Masks{i}.NL)
        continue;
    end
    
    disp(Dsdt(i).name);
    %load sdt file
    sdt = bh_readsetup([path Dsdt(i).name]);
    AcqRng = GetPhotonCollectTRange([path Dsdt(i).name]);
    timestp(i) = AcqRng(1);
    
    % Get indexes
    dashes = strfind(Dsdt(i).name,'_');
    
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
    filename = Dsdt(i).name;
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
    
    %time axis
    decay = zeros(size(flim,1),1);
    time = (1:length(decay))'*dt;
    
    fn = fieldnames(Masks{i}(1));
    if isempty(eval(['Masks{i}(1).' fn{1}]))
        continue;
    end
    
    % Masks has NADH-FAD structure, but there is only one mask for each sdt
    % in this workflow. Both 'NL' and 'FL' are the same
    selected_pixel=Masks{i}(1).NL;
    
    reshapedselectedpix = repmat(reshape(uint16(selected_pixel),[1,size(selected_pixel)]),[length(decay),1,1]);
    decay = decay + sum(sum(reshapedselectedpix.*flim,2),3);
    
    % Calculate irradiance and irradiance std for averages
    % later
    ind = find(selected_pixel);
    IntIm = squeeze(sum(flim,1));
    IntVals = IntIm(ind);
    photons = sum(IntVals);
    irr = mean(IntVals)/numscans;
    irr_std = std(IntVals)/numscans;
    if ProfBool
        IntImSc = IntIm./IllProfCal.*mean(IllProfCal(:));
        IntValsSc = IntImSc(ind);
        irrSc = mean(IntValsSc)/numscans;
        irrSc_std = std(IntValsSc)/numscans;
    else
        irrSc = -1;
        irrSc_std = -1;
    end
    
    %update total photon counts
    totphot = num2str(sum(decay));
    
    %% Save decays in the same format that is output by 'FLIMDataAnalysisGUI_ver5_3'
    decay_structs{i,1}.decay = decay;
    decay_structs{i,1}.name = filename;
    decay_structs{i,1}.filename = filename;
    decay_structs{i,1}.image = image;
    decay_structs{i,1}.selected_pixel = selected_pixel;
    decay_structs{i,1}.time = time;
    decay_structs{i,1}.num_pixel = sum(sum(selected_pixel));
    decay_structs{i,1}.numscans = numscans;
    decay_structs{i,1}.photons = photons;
    decay_structs{i,1}.irr = irr;
    decay_structs{i,1}.irr_std = irr_std;
    decay_structs{i,1}.irrSc = irrSc;
    decay_structs{i,1}.irrSc_std = irrSc_std;
    decay_structs{i,1}.timestp = timestp(i);
    
    fclose('all');
end

% Final save
decay_struct = decay_structs;
save([path 'decays_SingleMasks.mat'],'decay_struct');