function [ml,MDpars,txt] = GroupMasterListPars(listpaths,parpaths)
% Function for grouping together master lists and corresponding MDpars from
% multiple different days. E.g. if you want to plot replicates together.
% INPUTS:
% -listpaths: cell of paths to the master lists (typically 'MasterList.xls').
% -parpaths: cell of paths to the corresponding MDpars files (typically 'MasterList_pars.mat').
% OUTPUTS: concatenated table and array (text table optional)

% clear all;
% load('C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\temp')
% listpaths{1} = 'MasterList.xls'; listpaths{2} = 'MasterList2.xls';
% parpaths{1} = 'MasterList_pars.mat'; parpaths{2} = 'MasterList_pars2.mat';

if ~iscell(listpaths) listpaths2{1} = listpaths; listpaths = listpaths2; clear('listpaths2'); end
if ~iscell(parpaths) parpaths2{1} = parpaths; parpaths = parpaths2; clear('parpaths2'); end

% % Data set 1
% [num,txtcat,mlcat]= xlsread(listpaths{1});
% % Remove 1st row of headers
% mlcat(1,:)=[]; txtcat(1,:)=[];
% load(parpaths{1}); MDparscat = MDpars;
% % Dims: [Param#, mean/std(1,2), ListRow, time point, channel, z-position, segment]

% Loop through remaining lists and concatenate them to the first, along the
% row dimension.
ml=[];
for i=1:size(listpaths,2)
    [num,txt,ml0]= xlsread(listpaths{i});
    % Remove 1st row of headers
    ml0(1,:)=[]; txt(1,:)=[];
    mlClms = size(ml,2); ml0Clms = size(ml0,2);
    if mlClms > ml0Clms 
        ml0 = [ml0 num2cell(nan(size(ml0,1),mlClms-ml0Clms))];
    elseif mlClms < ml0Clms 
        ml = [ml num2cell(nan(size(ml,1),ml0Clms-mlClms))];
    end
        
    ml = [ml; ml0];
    
    load(parpaths{i});
    MDs{i} = MDpars;
    dims(i,:) = size(MDpars);
end

% To concat MDpars, they need to have the same dims. So find max dims and
% pad smaller arrays with NaN's to make them all the same size for
% concatenating along 'List row' dim.
mxdims = max(dims,[],1);
for i = 1:length(MDs)
    dims = size(MDs{i});
    % Find dims difference to pad mats with NaNs
    dimsdiff = mxdims-dims;
    % But don't pad 'row' dimension. That's the one we're concatenating along
    dimsdiff(3) = 0;
    MDs{i} = padarray(MDs{i},dimsdiff,nan,'post');
end

MDpars = cat(3,MDs{:});
% Omit rows with a '1' in the 'Exclude' column
omt = [ml{:,4}]==1;
ml(omt,:)=[]; MDpars(:,:,omt,:,:,:,:,:,:)=[]; 
% txt(omt,:)=[]; % Not really using this
