function FLIM_CatT_mdpars(acqs,MskCorrFile)
% Simple function for concatenative separate acquisition mdpars matrices
% along the time dimension, then replotting.
% Written for poison pipetting experiments, where we took baseline, data,
% then stopped the acq, pipetted poison, then started the 'part B'
% acquisition, but it's really one time course.
% NOTE: assume masks are numbered the same.
% INPUTS:
% -acqs: a cell of path names to acquisition folders to be concatenated
% -MskCorrFile: path to an xls spreadsheet containing mask number
%   correspondences. Useful if masks switch between acqs, e.g eggs move
%   E.g. column 1 has mask nums in its rows for first acq, then column 2
%   has mask nums for 2nd acq that correspond. Blank rows are filled with
%   nans.

% clear all;
% acqs{1} = 'Z:\Lab\Marta\Cumulus_cell_New_Scope\2019-07-31_Test45\s6_a1a\';
% acqs{2} = 'Z:\Lab\Marta\Cumulus_cell_New_Scope\2019-07-31_Test45\s6_a1b_Arsenite\';
% acqs{3} = 'Z:\Lab\Marta\Flow Chamber Tests\2019-05-07 Rotenone w oxamate pipette\2019-05-15 Third repeat\s1_a1c\';
% MskCorrFile = 'Z:\Lab\Marta\Flow Chamber Tests\2019-05-07 Rotenone w oxamate pipette\2019-05-15 Third repeat\s1_a1_msk_corr.xls';

% if ~iscell(acqs) acqs{1} = acqs; end
if acqs{1}(end)~='\'; acqs{1} = [acqs{1} '\']; end
startupTim

% Check if only 'JointMasks' were produced, then assume MaskLab = 'JointMasks'
Dmsk = dir([acqs{1} 'Masks*.mat']); Djnt = dir([acqs{1} 'JointMasks*.mat']);
MaskLab = '';
if isempty(Dmsk)&~isempty(Djnt) MaskLab = 'JointMasks'; end

load([acqs{1} 'multiD_pars_' MaskLab '.mat']);
load([acqs{1} 'multiD_indices.mat']); nameindscat = nameinds;
load([acqs{1} 'BG_vals.mat']);
if exist('MskCorrFile')&MskCorrFile~=-1 mskcorr = xlsread(MskCorrFile); end

% Find cal_files folder
if exist([acqs{1} 'cal_files'])==7 % if cal_files is in acqs{1}
    calpath = [acqs{1} 'cal_files\'];
elseif exist([UpOneDir(acqs{1}) 'cal_files'])==7 % if cal_files in daypath
    calpath  = [UpOneDir(acqs{1}) 'cal_files\'];
else
    error('Cannot locate cal_files folder. Place in daypath or acqpath, please');
end

for a = 1:size(acqs,2)
    if acqs{a}(end)~='\'; acqs{a} = [acqs{a} '\']; end
    load([acqs{a} 'multiD_pars_' MaskLab '.mat']);
    load([acqs{a} 'multiD_indices.mat']);
    load([acqs{a} 'BG_vals.mat']);
    dims(a,:) = size(mdpars);
    mds{a} = mdpars;
    
    nameindss{a} = nameinds;
    BGvalss{a} = BGvals;
end

% To use this function, positions really have to coordinate. If sub-acqs
% have different numbers of positions, error. Go back and use
% 'MultiD_exclude_frames.m' to remove positions that aren't present in ALL
% sub-acqs (and rename the Position folders to coordinate, if necessary),
% then re-run this cat function
if length(unique(dims(:,4)))>1
    error('Different numbers of positions detected. Sub-acquisitions must have corresponding positions.')
end

% To concat mdpars, they need to have the same dims. So find max dims and
% pad smaller arrays with NaN's to make them all the same size for
% concatenating along 'time' dim.
mxdims = max(dims,[],1);
for i = 1:length(mds)
    dims = size(mds{i});
    % Find dims difference to pad mats with NaNs
    dimsdiff = mxdims-dims;
    % But don't pad 'row' dimension. That's the one we're concatenating along
    dimsdiff(3) = 0;
    mds{i} = padarray(mds{i},dimsdiff,nan,'post');
end

% Make directories for copied tiffs
for p = 1:mxdims(4)
    [a,b] = mkdir([acqs{1}(1:end-1) '_cat\sorted_sdts\IntTiffs_Pos' num2str(p) '\']);
    % D{1} = dir([acqs{1} '\sorted_sdts\IntTiffs_*\*.tif']);
    [a,b] = mkdir([acqs{1}(1:end-1) '_cat\sorted_sdts\SHGTiffs_Pos' num2str(p) '\']);
    % DS{1} = dir([acqs{1} '\sorted_sdts\SHGTiffs*\*.tif']);
    [a,b] = mkdir([acqs{1}(1:end-1) '_cat\sorted_sdts\ROIsCheck_Pos' num2str(p) '\']);
    % DR{1} = dir([acqs{1} '\sorted_sdts\SHGTiffs*\*.tif']);
    dims(1,:) = size(mdpars);
end

% Mask resorting
if exist('mskcorr')
%     rm = []; % Flag masks that went missing for removal... mmm actually
%     not necessary. Just exclude later with plotting if we want.
    for a = 2:size(acqs,2)
        if a>1
%             rm = [rm find(isnan(mskcorr(:,a)))];
            mdsrt = nan(size(mds{a}));
            for i = 1:mxdims(7) % Go down rows and fill in sorted array if index is filled in. 
                if ~isnan(mskcorr(i,1))&~isnan(mskcorr(i,a))
                    mdsrt(:,:,:,:,:,:,mskcorr(i,1),:)=mds{a}(:,:,:,:,:,:,mskcorr(i,a),:);
                end
            end
%             [(1:mxdims(7))' squeeze(mds{a}(6,1,1,1,1,1,:,1)) squeeze(mdsrt(6,1,1,1,1,1,:,1))] % test
            mds{a} = mdsrt;
        end
    end
end
mdparscat = cat(3,mds{:} );
nameindscat = cat(1,nameindss{:});
BGvalscat = cat(2,BGvalss{:});
% BGnum_pixelscat = cat(2,BGnum_pixelss{:}); BGstdscat = cat(2,BGstdss{:}); BGtimestpscat = cat(2,BGtimestpss{:});

for a = 1:size(acqs,2)
for p = 1:mxdims(4) % Note, already checked above for sub-acqs with different nums of positions
    % Int tiffs directories
    D{a,p} = dir([acqs{a} '\sorted_sdts\IntTiffs_Pos' num2str(p) '*\*.tif']);
    DS{a,p} = dir([acqs{a} '\sorted_sdts\SHGTiffs_Pos' num2str(p) '*\*.tif']);
    DR{a,p} = dir([acqs{a} '\sorted_sdts\ROIsCheck_Pos' num2str(p) '*\*.png']);
    
    for i = 1:length(D{a,p})
        dest = [acqs{1}(1:end-1) '_cat\sorted_sdts\IntTiffs_Pos' num2str(p) '\Acq' num2str(a) '_' D{a,p}(i).name];
        [a0,b] = copyfile([D{a,p}(i).folder '\' D{a,p}(i).name],dest);
    end
    for i = 1:length(DS{a,p})
        dest = [acqs{1}(1:end-1) '_cat\sorted_sdts\SHGTiffs_Pos' num2str(p) '\Acq' num2str(a) '_' D{a,p}(i).name];
        [a0,b] = copyfile([DS{a,p}(i).folder '\' DS{a,p}(i).name],dest);
    end
    for i = 1:length(DR{a,p})
        dest = [acqs{1}(1:end-1) '_cat\sorted_sdts\ROIsCheck_Pos' num2str(p) '\Acq' num2str(a) '_' D{a,p}(i).name(1:end-3) 'png'];
        [a0,b] = copyfile([DR{a,p}(i).folder '\' DR{a,p}(i).name],dest);
    end
end
end
% Rename final files and save to output dir
mdpars = mdparscat;
nameinds = nameindscat;
BGvals = BGvalscat;
save([acqs{1}(1:end-1) '_cat\' 'multiD_pars_' MaskLab '.mat'],'mdpars');
save([acqs{1}(1:end-1) '_cat\' 'multiD_indices.mat'],'nameinds');
save([acqs{1}(1:end-1) '_cat\' 'BG_vals.mat'],'BGvals');
copyfile(calpath,[acqs{1}(1:end-1) '_cat\cal_files']);

% Also copy the fits*.mat files from the 1st acquisition, so that labels
% and some metadata can be used by subsequent programs like
% FLIM_master_list_and_MDpars.m
% NOTE: assume each sub-folder has the same number of fits, and that mask
% correspondence has been properly accounted for with 'MskCorrFile'
Dfits = dir([acqs{1} '*fits*.mat']);
for i = 1:length(Dfits)
    copyfile([Dfits(i).folder '\' Dfits(i).name],[acqs{1}(1:end-1) '_cat\' Dfits(i).name]);
end

LaserOPplot(acqs,[acqs{1}(1:end-1) '_cat\LaserOPs'])

ts = unique(squeeze(mdparscat(7,1,:,:,:,:,:,:,:,:))); ts(isnan(ts))=[];
if length(ts)>1 FLIMAcqParamPlot([acqs{1}(1:end-1) '_cat']); end

