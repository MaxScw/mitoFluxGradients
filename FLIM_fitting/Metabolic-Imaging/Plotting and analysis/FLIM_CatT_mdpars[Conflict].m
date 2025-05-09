function FLIM_CatT_mdpars(acqs)
% Simple function for concatenative separate acquisition mdpars matrices
% along the time dimension, then replotting.
% Written for poison pipetting experiments, where we took baseline, data,
% then stopped the acq, pipetted poison, then started the 'part B'
% acquisition, but it's really one time course.
% NOTE: assume masks are numbered the same.

% acqs{1} = 'C:\Dropbox\temp\s1_a1\';
% acqs{2} = 'C:\Dropbox\temp\s1_a1b\';
% acqs{3} = 'Z:\Lab\Marta\Flow Chamber Tests\2019-04-24 Oxamate FCCP_pipette\s1_a1c\';
% acqs{4} = 'Z:\Lab\Marta\Flow Chamber Tests\2019-04-24 Oxamate FCCP_pipette\s1_a1d\';

% if ~iscell(acqs) acqs{1} = acqs; end
if acqs{1}(end)~='\'; acqs{1} = [acqs{1} '\']; end

load([acqs{1} 'multiD_pars.mat']); 
load([acqs{1} 'multiD_indices.mat']); nameindscat = nameinds;
load([acqs{1} 'BG_masks_ints.mat']); BGirrscat = BGirrs; % BGmaskscat = BGmasks; BGnum_pixelscat = BGnum_pixels; BGstdscat = BGstds; BGtimestpscat = BGtimestps; 
[a,b] = mkdir([acqs{1}(1:end-1) '_cat\IntTiffs\']);
D{1} = dir([acqs{1} '\sorted_sdts\IntTiffs_*\*.tif']);
[a,b] = mkdir([acqs{1}(1:end-1) '_cat\SHGTiffs\']);
DS{1} = dir([acqs{1} '\sorted_sdts\SHGTiffs*\*.tif']);
dims(1,:) = size(mdpars);

for a = 1:size(acqs,2)
    if acqs{a}(end)~='\'; acqs{a} = [acqs{a} '\']; end
    load([acqs{a} 'multiD_pars.mat']);
    load([acqs{a} 'multiD_indices.mat']);
    load([acqs{a} 'BG_masks_ints.mat']); 
    dims(a,:) = size(mdpars);
    mds{a} = mdpars;
    
    nameindss{a} = nameinds;
    BGirrss{a} = BGirrs;
%     nameindscat = [nameindscat;nameinds];
    

    
    % Int TiffsF
    D{a} = dir([acqs{a} '\sorted_sdts\IntTiffs_*\*.tif']);
    DS{a} = dir([acqs{a} '\sorted_sdts\SHGTiffs*\*.tif']);
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
mdparscat = cat(3,mds{:} );
nameindscat = cat(1,nameindss{:});
BGirrscat = cat(2,BGirrss{:});


frind = 1; % Index for copying IntTiffs and SHGTiffs
for a = 1:size(D,2)
        for i = 1:length(D{a})
            dest = [acqs{1}(1:end-1) '_cat\IntTiffs\fr' num2str(frind,'%05i') '.tif'];
            [a0,b] = copyfile([D{a}(i).folder '\' D{a}(i).name],dest);
            dest = [acqs{1}(1:end-1) '_cat\SHGTiffs\fr' num2str(frind,'%05i') '.tif'];
            [a0,b] = copyfile([DS{a}(i).folder '\' DS{a}(i).name],dest);
            frind = frind + 1;
        end
end
mdpars = mdparscat;
nameinds = nameindscat;
BGirrs = BGirrscat;
save([acqs{1}(1:end-1) '_cat\' 'multiD_pars.mat'],'mdpars');
save([acqs{1}(1:end-1) '_cat\' 'multiD_indices.mat'],'nameinds');
save([acqs{1}(1:end-1) '_cat\' 'BG_masks_ints.mat'],'BGirrs');

FLIMAcqParamPlot([acqs{1}(1:end-1) '_cat'])
