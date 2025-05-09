function MultiD_sdt_sort_uMan_OBS(path)
% Using uManager, multi-D acquisitions can be taken with any of the
% following dimensions: [Pos,t,channel,z]
% When saving B&H sdt's from the ttl trigger from the z piezo, sdt files
% are saved in one big series. You also save the uManager demo cam images
% in the same folder as the multi-D acquisition (use 4x4 binning to reduce
% image size).
% This program uses the uManager structure and image time stamps to sort
% and rename the sdt's for further analysis.

% clear all; path = 'C:\Users\Tim\Documents\Academic - Research\Data\2017-06-29 Labview test acq\a3\';

% Version:
% 2017-07-01: Added case for handling acquisitions taken with Labview.
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
D = dir([path 'uMan_*']);
uManpath = [path D(end).name '\'];
Dsdt = dir([path '*.sdt']);
sdtL = size(Dsdt,1);
[a,b] = mkdir(sdtpath);
DMan = dir(uManpath); DMan(1:2)=[]; DMan(~[DMan.isdir])=[];

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


for posnum = 1:max([length(DMan),1])
    % Assume at least 1 position. If no pos specified, assume only 1 pos and
    % tiffs are in uManpath.
    if isempty(DMan)
        DManfiles{1} = dir([uManpath '\*.tif']);
        DManmeta{1} = dir([uManpath '\metadata.txt']);
        sdtsortedpaths{posnum}=[sdtpath 'Pos0\'];
        uManPos = 0;
    else
        % Get subdirectories
        DManfiles{posnum} = dir([uManpath DMan(posnum).name '\*.tif']);
        DManmeta{posnum} = dir([uManpath DMan(posnum).name '\metadata.txt']);
        sdtsortedpaths{posnum}=[sdtpath DMan(posnum).name '\'];
    end
    % Was this acquisition a uMan 'zscan'? If so, there will be at least
    % one '001' in the last index in the uMan file. Then use the last index
    % for the z instead of 'z0,z1', etc.
    if isempty(strfind([DManfiles{1}.name],'_z0_on_001'))
        zScanBool = 0;
    else
        zScanBool = 1;
    end
    
    % Get rid of 'waiting' files used for laser wavelength switching
    remove = [];
    for k = 1:length(DManfiles{posnum})
        if ~isempty(strfind(DManfiles{posnum}(k).name,'wait'))|~isempty(strfind(DManfiles{posnum}(k).name,'_off_'))|~isempty(strfind(DManfiles{posnum}(k).name,'trg'))|~isempty(strfind(DManfiles{posnum}(k).name,'temp'))|~isempty(strfind(DManfiles{posnum}(k).name,'tbuff'))|~isempty(strfind(DManfiles{posnum}(k).name,'trans'))
            remove = [remove k];
        end
%         if ~isempty(strfind(DManfiles{posnum}(k).name,'wait'))
%             remove = [remove k];
%         end
    end
    DManfiles{posnum}(remove)=[];
    
    % Find earliest non-waiting uMan file
    AllManFiles = cell2mat([DManfiles(:)]);
%     earliestuManind = find([AllManFiles(:).datenum] == min([DManfiles{1}(:).datenum]));
    tstampoffset = min([AllManFiles(:).datenum])-AcqTRng1(2); % stop time of photon coll.
    [a,b] = mkdir(sdtsortedpaths{posnum});
end

for posnum = 1:max([length(DMan),1])
    prevSdtInd = -1;
    for i = 1:min([size(Dsdt,1),length(DManfiles{posnum})])
        clear tstamps_uMan
        % These uMan timestamps only refer to the save time of the uMan
        % file. Since these files are tiny, they save instantly, so this
        % time should be very close to when the stop time of photon
        % collection.
        
%         if i ==160
%             disp('')
%         end
        
        file = DManfiles{posnum}(i).name;
        % DAQ trigger used, channels will have '_z0_'
        dashes = strfind(file,'_');
        t = str2num(file(dashes(1)+1:dashes(2)-1));
        chan = file(dashes(2)+1:dashes(3)-1);
        if zScanBool
            z = str2num(file(dashes(5)+1:end-4));
        else
            z = str2num(file(dashes(3)+2:dashes(4)-1));
        end
        tstamps_uMan(i) = DManfiles{posnum}(i).datenum;
        % Using time stamps find the corresponding sdt file, which is
        % the file closest in time to current uMan file.
        % sdt should save shortly before uMan file does.
        
        diffs = tstamps_sdt-tstamps_uMan(i)+tstampoffset;
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
        dest = [sdtsortedpaths{posnum} 'sdt_' num2str(t,'%09g') '_' chan '_' num2str(z,'%03g') '.sdt'];
        Dtemp = dir(dest);
        if isempty(Dtemp) % Copy if file hasn't already been copied.
            copyfile(orig,dest);
        end

        % Add the frame to the name index. Include the uMan timestamp
        % in case of duplicates, but delete that column at the end
        % Use Columns 7 and 8 (timestamps and nscans) as the flag column.
        % This keeps sdt, z, t, and tiff correspondences, but timestamp
        % won't be needed.
        if ~isempty(DMan) 
            uManPos = DMan(posnum).name(4:end);
%             strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
%             uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
        end
        % Get number of scans to later make scaled intensity matrix
        sdt = bh_readsetup(orig);
        meas = bh_getmeasdesc(sdt, 1);
        numscans = meas.hist_fida_points;
        
        % For uMan positions with arbitrary string labels (not 'Pos_')
%         if isempty(uManPos) uManPos = DMan(posnum).name; end
        nameinds = [nameinds; {SdtMatchingInd,t,uManPos,chan,z,-1,AcqTRng{SdtMatchingInd}(1),AcqTRng{SdtMatchingInd}(2),numscans}];
         
        if ~isempty(dups)
            diff1 = abs(tstamps_sdt(SdtMatchingInd)-nameinds{dups,8});
            diff2 = abs(tstamps_sdt(SdtMatchingInd)-tstamps_uMan(i));
            if diff1<diff2
                nameinds(end,[7 8 9])=num2cell(-1); % Flag the row just added, ie, current frame
            else
                nameinds(dups,[7 8 9])=num2cell(-1); % Flag other duplicate
            end
        end
    end
end

nameinds = sortrows(nameinds,1);
save([path 'multiD_indices.mat'],'nameinds');


