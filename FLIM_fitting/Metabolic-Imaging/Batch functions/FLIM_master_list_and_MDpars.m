function FLIM_master_list_and_MDpars(daypath,ExcludeFolders,Label)
% Once everything has been calculated for a big batch, use this function to
% make a master list of all fits from a single day (root dir). This list
% can be manually edited to set plot groups and exclusions.
% Function also re-orders mdpars to simplify the dimensionality. Pos# and
% mask# get collapsed into a single dimension that is a linear list of
% fits the same correspondence to the master list. New array is 'MDpars'.
% -ExcludeFolders: Sometimes you don't want to analyze all folders with
%    fits in them. ExcludeFolders is a cell of folder names to exclude.
% -Label: if you want to create a custom label for the output averages file

% % Testing
% clear all;
% daypath = 'Z:\Lab\Tim\Boston IVF\Discarded_study\Marta\2019-10-11_TestAcqui\';
% % ExcludeFolders = {'Pre201809Revision'};

if daypath(end)~='\' daypath = [daypath '\']; end;
if ~exist('ExcludeFolders')|~iscell(ExcludeFolders) clear ExcludeFolders; end % ExcludeFolders must be a cell of paths, even if there is only one element
if ~exist('Label')|Label==-1 clear Label; end
close all;

sams = {};
Dfits = [];
fnames = {};
pos = {};
mask = {};
labels = {};

% Do a quick for loop. Go folder by folder (presumably sample by sample),
% get positions and sort by position. Build a 'Dfits' structure in that
% order.
Df = dir(daypath); Df(1:2)=[]; Df(~[Df.isdir])=[];
if exist('ExcludeFolders')
    Exc=[];
    for i = 1:length(Df)
        for j = 1:length(ExcludeFolders)
            if strcmp(Df(i).name,ExcludeFolders{j})
                Exc = [Exc i];
            end
        end
    end
    Df(Exc)=[];
end

for i = 1:length(Df)
    Dfits_tmp = dir([daypath Df(i).name '\*fits*.mat']);
    [a,b] = natsort({Dfits_tmp.name});
    Dfits_tmp = Dfits_tmp(b); % Natural-sort by positions
    Dmdpars = dir([daypath Df(i).name '\*multiD_pars*.mat']);
    Dsorted_sdts = dir([daypath Df(i).name '\sorted_sdts']);
    % Skip folders that either don't have multiD_pars.mat or a sorted_sdts
    % folder or fits
    if isempty(Dmdpars)|isempty(Dsorted_sdts)|isempty(Dfits_tmp)
        continue;
    end
    clear sams_tmp fnames_tmp pos_tmp mask_tmp labels_tmp
    for j = 1:length(Dfits_tmp)
        sams_tmp{j} = Df(i).name;
        fnames_tmp{j} = Dfits_tmp(j).name;
        dashes = strfind(fnames_tmp{j},'_');
        posind = strfind(fnames_tmp{j},'Pos');
        %         maskind = strfind(fnames_tmp{j},'mask');
        MasLab = fnames_tmp{j}(dashes(2)+1:strfind(fnames_tmp{j},'.mat')-1);
        
        if ~isempty(posind)
            pos_tmp{j} = (fnames_tmp{j}(posind+3:dashes(2)-1));
            %             mask_tmp{j} = fnames_tmp{j}(maskind+4:dashes(end)-1);
            mask_tmp{j} = MasLab;
            if isempty(strfind(fnames_tmp{j},'thr'))
                labels_tmp{j} = [sams_tmp{j} 'P' pos_tmp{j} mask_tmp{j}];
            else
                labels_tmp{j} = [sams_tmp{j} 'P' pos_tmp{j} 'thr'];
            end
        else
            pos_tmp{j} = -1;
            mask_tmp{j} = -1;
            labels_tmp{j} = [fnames_tmp{j}];
        end
    end
    
    sams = [sams sams_tmp];
    Dfits = [Dfits; Dfits_tmp];
    fnames = [fnames fnames_tmp];
    pos = [pos pos_tmp];
    mask = [mask mask_tmp];
    labels = [labels labels_tmp];
end
usams = unique(sams);
% Load mdpars and get indices to use the right mdpars for each acquisition

% Calculate dimensions for MDpars by looping through samples
szs = [];
for i = 1:length(usams)
    % May be more than one multiD_pars. Load any one of them to get the
    % dims for this sample.

    % Get dimensions of MDpars to pre-initialize to NaNs. Use mdpars from
    % acq folders
    % MDpars dims: [1=Param#, 2=mean/std(1,2), 3=ListRow, 4=time point, 5=channel, 6=z-position, 7=segment]
    % mdpars: [1=Param#, 2=mean/stderr(1,2), 3=time point, 4=position, 5=channel, 6=z-position, 7=embryo number, 8=segment]
    Dmdpars = dir([daypath usams{i} '\multiD_pars*.mat']);
    for mdi = 1:length(Dmdpars)
        load([Dmdpars(mdi).folder '\' Dmdpars(mdi).name]);
        szs = [szs; size(mdpars)];
        clear mdpars;
    end
    % Find bin for this sample
    Dec1InSam = find(strcmp(sams,usams{i})); Dec1InSam = Dec1InSam(1);
    load([daypath '\' usams{i} '\' Dfits(Dec1InSam).name]) % loads decays_fits_struct
    tmp = decays_fits_struct{~cellfun('isempty',decays_fits_struct)}; Bin = tmp.Bin;
    if szs(i,6)==Bin szs(i,6) = 1; end % If only one data point per time point, put in 1st element of z-dim
    Segdim(i) = size(tmp.decay,2);
end

% Initialize MDpars to NaNs.
szs = max(szs,[],1); % Get max dim sizes
% NOTE: channels should always have 3 elements now: 1=NADH, 2=FAD, 3=UserChan
MDpars = nan(max(szs(1)),max(szs(2)),size(Dfits,1),max(szs(3)),3,max(szs(6)),max(Segdim));

% For spreadsheets (one for means, one for std's)
spr{1,1} = 'Sample'; spr_std{1,1} = 'Sample';
% In the spreadsheet, use the following fields for manual labeling
spr{1,2,1} = 'Data Set';
spr{1,3,1} = 'Set Name';
spr{1,4,1} = 'Exclude(enter 1)';
spr{1,5,1} = 'Point Text';

for i =1:size(Dfits,1)
    disp([num2str(i) '/' num2str(size(Dfits,1))]);
    % Load fits and get metadata
    fname = Dfits(i).name;
    load([daypath sams{i} '\' fname]) % loads decays_fits_struct
    try     load([daypath sams{i} '\' 'multiD_indices.mat']); catch     load([daypath sams{i} '\' 'name_indexes.mat']); end
    PosNames = unique(nameinds(:,3));
    dashes = strfind(fname,'_');
    Pstr = fname(dashes(1)+4:dashes(2)-1);
    Pind = find(strcmp(PosNames,Pstr));
    Zdim = length(unique([nameinds{:,5}]));
    if length(dashes)>2
        Mind = str2num(fname(dashes(2)+5:dashes(3)-1)); % Not really used anymore
    else
        Mind = str2num(fname(dashes(2)+5:end-4));
        MaskLab = fname(dashes(2)+1:end-4);
    end
    if isempty(Mind) Mind = 1; end % Assume single mask if 'mask' not present
    
    % Enter name into list
    spr{i+1,1} = labels{i};
    
    % Select appropriate mdpars array
    if strcmp(MaskLab(1:4),'mask')
        load([daypath sams{i} '\multiD_pars.mat'])
    else
        load([daypath sams{i} '\multiD_pars_' MaskLab '.mat'])
    end
    
    if strcmp(fname,'fits_Pos0_mask8.mat')
        1;
    end
    
    % Load mdpars values for this decay (ie, pos and mask) into the master
    % MDpars matrix. Do so for each row, or pos-mask combo
    mdparsrw = mdpars(:,:,:,Pind,:,:,Mind,:); % mdpars for this row, ie sample
    sz = size(mdparsrw);
    % Use reshape to remove the Mind singleton index, then permute to match
    mdparsrw = permute(reshape(mdparsrw,[sz(1) sz(2) sz(3) sz(4) sz(5) sz(6) sz(8) ]),[1 2 4 3 5 6 7]);
    sz = size(mdparsrw);
    
    % load MDpars dimensions (for row 'i')
%     size(MDpars)
%     if i==10
%         1;
%     end
    MDpars(:,:,i,1:sz(4),:,1:sz(6),1:sz(7)) = mdparsrw;
    
    % Get embryo label for plotting
    dashes = strfind(Dfits(i).name,'_');
    PosInd = strfind(Dfits(i).name,'Pos');
    MaskInd = strfind(Dfits(i).name,'mask');
    PosNum = Dfits(i).name(PosInd+3:dashes(2)-1);
    if length(dashes)>2 % E.g. '..._fxshft.mat'
        MaskNum = Dfits(i).name(MaskInd+4:dashes(3)-1);
    else
        MaskNum = Dfits(i).name(MaskInd+4:end-4);
    end
    
    
    % Changepoints: search and load into spr, if present
    Dcp = dir([Dfits(i).folder '\embryo_data_changepoints.mat']);
    if ~isempty(Dcp)
        spr{1,6,1} = '1st Div';
        spr{1,7,1} = '2nd Div';
        spr{1,8,1} = 'Blastocoel';
        load([Dcp.folder '\' Dcp.name])
        Pcpind = strfind(chngpnts_accum(:,1),['Pos' PosNum]);
        Pcpind = ~cellfun('isempty',Pcpind);
        Mcpind = strfind(chngpnts_accum(:,1),['mask' MaskNum]);
        Mcpind = ~cellfun('isempty',Mcpind);
        cpind = find(Pcpind&Mcpind);
        cps = chngpnts_accum{cpind,2};
        if length(cpind)>1
            error('Unique labeling of changepoints violated.')
        end
        % Note: the frames in chngpnts_accum are the time points, and
        % that's actually correct. FLIM_time_plots.m uses those same units
        % (ie, 'time point', not 'frame num'
        spr{i+1,6} = cps(1);matlab individual license
        spr{i+1,7} = cps(2);
        spr{i+1,8} = cps(3);
    end
    
    dashes = strfind(sams{i},'_');
    if strcmp(sams{i}(1),'s')
        if isempty(dashes)
            spr{i+1,5} = [sams{i} 'P' PosNum 'm' MaskNum];
        else
            spr{i+1,5} = [sams{i}(1:dashes(1)-1) 'P' PosNum 'm' MaskNum];
        end
    else
        spr{i+1,5} = ['P' PosNum 'm' MaskNum];
    end
end

if exist('Label')
    Label = ['_' Label];
else
    Label = '';
end

% If previous MasterList was written, check if rows, i.e. samples are
% same? If so, the samples haven't changed, so there is no need to
% overwrite the MasterList, deleting previous annotations.
if exist([daypath '\MasterList' Label '.xls'])==2 % File exists
    [num,txt,raw] = xlsread([daypath '\MasterList' Label '.xls']);
    if length(spr(:,1))~=length(raw(:,1)) | isempty(find(strcmp(spr(:,1),raw(:,1)))) % rows have changed
        % Make a copy of previous MasterList so you don't accidentally
        % lose annotations.
        copyfile([daypath '\MasterList' Label '.xls'],[daypath '\MasterList' Label datestr(datetime,'yyyymmddhhss') '.xls'])
        delete([daypath '\MasterList' Label '.xls']);
        xlswrite([daypath '\MasterList' Label '.xls'],spr);
    end
else
    xlswrite([daypath '\MasterList' Label '.xls'],spr);
end
% Re-write MDpars no matter what, assuming we're running this because data
% has been updated somehow.
save([daypath 'MasterList_pars' Label '.mat'],'MDpars');
