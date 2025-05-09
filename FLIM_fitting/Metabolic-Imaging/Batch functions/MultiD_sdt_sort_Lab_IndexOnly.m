function MultiD_sdt_sort_Lab_IndexOnly(path)
% Using Labview and BH, multi-D acquisitions can be taken with any of the
% following dimensions: [Pos,t,channel,z]
% When saving B&H sdt's from the ttl trigger from the z piezo, sdt files
% are saved in one big series. Labview saves metadata in a spreadsheet
% in the same folder as the multi-D acquisition.

% This program uses the metadata and time stamps to sort and rename the
% sdt's for further analysis.

% clear all; path = 'Z:\Lab\Tim\Misc Tests\2017-06-29 Labview test acq\a6\';

% Version:
% 2017-07-01: Adapted for handling acquisitions taken with Labview.
%  Reads a spreadsheet text file for datestamps instead of uMan dummy tiffs
% 2016-02-05: updated identification of collection time and matching with
%    uMan files
% -Just switched to using uMan to trigger B&H with a DAQ card
% because z piezo was inconsistent. Filenames are different, using a
% different channel for each z position. Change the reading scheme in this
% function, but keep the previous sdt-writing scheme so I don't have to
% change the subsequent programs.
% read:   'img_000000000_FAD_z0_on_000.tif'
% write:  'sdt_000000000_FAD_000.sdt'
% 2016-07-19: Sometimes I take z-scans with uManager. They use the last
% index in the filename for z. Allow for these acq types with an 'if'

if path(end)~='\'; path = [path '\']; end;

slashes = strfind(path,'\');
Run = path(slashes (end-1)+1:end-1);
sdtpath = [path 'sorted_sdts\'];
AcqSdtCorr = ReadTxtSprd2Cell([path 'AcqSdtCorr.txt']);
Dsdt = dir([path '*.sdt']);
sdtL = size(Dsdt,1);
[a,b] = mkdir(sdtpath);

% Make an list, 'nameinds' that keeps track of the exact file correspondence
% [sdt# , uManfile-Pos , uManfile-t , uManfile-chan , uManfile-z, AcqSt, AcqEnd, NumScans]
nameinds = {};

% AcqTRng is the start and stop times of the photon collection in
% datenum units (absolute time)
AcqTRng1 = GetPhotonCollectTRange([path Dsdt(1).name]);
for i = 1:length(Dsdt)
    AcqTRng{i} = GetPhotonCollectTRange([path Dsdt(i).name]);
    tstamps_sdt(i) = AcqTRng{i}(2);
end

Positions = unique(cell2mat(AcqSdtCorr(:,2)));

for posnum = 1:max([length(Positions),1])
    % Find earliest trigger datenum, calculate offset with save time of sdt
    tstampoffset = min(datenum(AcqSdtCorr(1,6)))-AcqTRng1(2); % stop time of photon coll.
    [a,b] = mkdir([sdtpath 'Pos' num2str(Positions(posnum)) '\']);
end

for i = 1:min([size(Dsdt,1),size(AcqSdtCorr,1)])
    clear tstamps_Lab
    
    % Get metadata info for i-th trigger
    t = AcqSdtCorr{i,3};
    chan = AcqSdtCorr{i,4};
    z = AcqSdtCorr{i,5};
    tstamps_Lab(i) = datenum(AcqSdtCorr{i,6});
    Pos = num2str(AcqSdtCorr{i,2});
    
    % Using time stamps find the corresponding sdt file, which is
    % the file closest in time to trigger.
    diffs = tstamps_sdt-tstamps_Lab(i)+tstampoffset;
    negdiffs = find(diffs<0);
    smallestdiffs = find(abs(diffs)==min(abs(diffs)));
    SdtMatchingInd = smallestdiffs;
    
    % Every once and while, you get two StdMatchingInd's, so take
    % the first one (earlier)
    SdtMatchingInd = SdtMatchingInd(1);
    % Sometimes delays happen because B&H takes too long to save. Sometimes
    % skip a frame. If that happens, two uMan files will get
    % matched with one sdt.
    % Solution: find the uMan that the sdt is closer to, flag the other
    % frame to be exluded from the rest of analysis.
    if ~isempty(nameinds)
        dups = find([nameinds{:,1}]==SdtMatchingInd&[nameinds{:,8}]>0);
    else
        dups = [];
    end
    orig = [path Dsdt(SdtMatchingInd).name];
    dest = [sdtpath 'Pos' Pos '\sdt_' num2str(t,'%09g') '_' chan '_' num2str(z,'%03g') '.sdt'];
    Dtemp = dir(dest);
    if isempty(Dtemp) % Copy if file hasn't already been copied.
        copyfile(orig,dest);
    end
    
    % Add the frame to the name index. Include the uMan timestamp
    % in case of duplicates, but delete that column at the end
    % Use Columns 7 and 8 (timestamps and nscans) as the flag column.
    % This keeps sdt, z, t, and tiff correspondences, but timestamp
    % won't be needed.
    
    %             strnums = sscanf(LabPos ,'%g'); %Find the numbers in the name
    %             LabPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    
    % Get number of scans to later make scaled intensity matrix
    sdt = bh_readsetup(orig);
    meas = bh_getmeasdesc(sdt, 1);
    numscans = meas.hist_fida_points;
    
    % For uMan positions with arbitrary string labels (not 'Pos_')
    %         if isempty(LabPos) LabPos = DLab(posnum).name; end
    nameinds = [nameinds; {SdtMatchingInd,Pos,t,chan,z,-1,AcqTRng{SdtMatchingInd}(1),AcqTRng{SdtMatchingInd}(2),numscans}];
    
    if ~isempty(dups)
        diff1 = abs(tstamps_sdt(SdtMatchingInd)-nameinds{dups,8});
        diff2 = abs(tstamps_sdt(SdtMatchingInd)-tstamps_Lab(i));
        if diff1<diff2
            nameinds(end,[7 8 9])=num2cell(-1); % Flag the row just added, ie, current frame
        else
            nameinds(dups,[7 8 9])=num2cell(-1); % Flag other duplicate
        end
    end
end

nameinds = sortrows(nameinds,1);
save([path 'name_indexes.mat'],'nameinds');


