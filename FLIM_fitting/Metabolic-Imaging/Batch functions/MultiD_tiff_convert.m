function MultiD_tiff_convert(acqpath, SHGBool, WellCropBool, bit16Bool, FLIMageBool)
% Take sdt files sorted by 'MultiD_sdt_sort', get integrated intensity
% images by summing along time domain, save them as tiffs in a subfolder.
% Inputs:
% -acqpath: acqpath to sample folder
% -SHGBool: '1' to convert SHG images in ch2 of sdt file
% -WellCropBool: '1' to look for circular edges of microwells and crop them out
% -bit16Bool: '1' to export 16bit tifs. Only necessary for very bright
%             samples (pixels get >255 photons)
% -FLIMageBool: '1' to export mean lifetime images with QuickFLIMage.m

% 1 tiff folder for each Pos folder that has all images in that folder in order of z, t
% If it is 1-chan data, just z slices of that one channel
% If it is 2-chan data (NADH and FAD), then concatenate the image pairs side-by-side
% acqpath - the acqpath to the main acquisition folder

% % TEST in script mode:
% clear all;
% acqpath = 'Z:\Lab\Marta\2016-11-17 ROS assay tests\2017-07-15 Batch 9\s1_a1B_Pos2Only_2FADfrs\';

if acqpath(end)~='\'; acqpath = [acqpath '\']; end;
if ~exist('SHGBool')|SHGBool==-1 SHGBool = 1; end
if ~exist('bit16Bool')|bit16Bool==-1 bit16Bool = 0; end
if ~exist('WellCropBool')|WellCropBool==-1 WellCropBool = 0; end
if ~exist('FLIMageBool')|FLIMageBool==-1 FLIMageBool = 0; end
slashes = strfind(acqpath,'\');
Run = acqpath(slashes (end-1)+1:end-1);
sdtpath = [acqpath 'sorted_sdts\'];
Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    if ~strcmp(Dpos(i).name(1:3),'Pos')
        remove = [remove i];
    end
end
Dpos(remove)=[];
block_1=1;
block_2=2;
% load nameinds. If frames were partially converted before, Stack frame
% correspondence in column 6 should still be there, so function will pick
% up where it left off.
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end

for posnum = 1:length(Dpos)
    Dnad = dir([sdtpath Dpos(posnum).name '\*NADH*.sdt']);
    Dfad = dir([sdtpath Dpos(posnum).name '\*FAD*.sdt']);
    Duser = dir([sdtpath Dpos(posnum).name '\*UserChan*.sdt']);
    if isempty(Dnad)&isempty(Dfad)&isempty(Duser) Duser = dir([sdtpath Dpos(posnum).name '\*Custom*.sdt']); end
    D = dir([sdtpath Dpos(posnum).name '*.sdt']);
    sdtL = max([length(Dnad) length(Dfad) length(Duser) length(D)]);
    tifpath = [sdtpath 'IntTiffs_' Dpos(posnum).name '_' Run '\'];
    FLIMagepath = [sdtpath 'FLIMageTiffs_' Dpos(posnum).name '_' Run '\'];
    SHGpath = [sdtpath 'SHGTiffs_' Dpos(posnum).name '_' Run '\'];
    [a,b] = mkdir(tifpath);
    if FLIMageBool [a,b] = mkdir(FLIMagepath); end
    [a,b] = mkdir(SHGpath);
    
    
    Dtiffs = dir([tifpath '*.tif']);
    PosNum = Dpos(posnum).name(4:end);
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    
    % See if a tiff conversion was started before. If so, determine which
    % tiffs have not been converted yet and convert those.
    % ALSO, refill in the correspondence in 'nameinds'. If it errored out
    % previously, nameinds may not have been saved.
    frs2do = 1:sdtL;
    inds=cell(sdtL,1);
    if ~isempty(Dtiffs)
        for k = 1:length(Dtiffs)
            frnum = str2num(Dtiffs(k).name(3:7));
            frs2do(frs2do==frnum)=[];
            
            if ~isempty(Dnad) % Only NADH images taken
                dashes = strfind(Dnad(frnum).name,'_');
                t = str2num(Dnad(frnum).name(dashes(1)+1:dashes(2)-1));
                z = str2num(Dnad(frnum).name(dashes(3)+1:end-4));
                %                 if t==5 & z ==1
                %                     disp('')
                %                 end
            elseif ~isempty(Dfad)
                dashes = strfind(Dfad(frnum).name,'_');
                t = str2num(Dfad(frnum).name(dashes(1)+1:dashes(2)-1));
                z = str2num(Dfad(frnum).name(dashes(3)+1:end-4));
            elseif ~isempty(Duser)
                dashes = strfind(Duser(frnum).name,'_');
                t = str2num(Duser(frnum).name(dashes(1)+1:dashes(2)-1));
                z = str2num(Duser(frnum).name(dashes(3)+1:end-4));
            else
                error('Channel error. No sdt''s detected.')
            end
            inds{frnum} = (PosInd)&([nameinds{:,2}]==t)&([nameinds{:,5}]==z);
        end
    end
    
    % IF DUAL CHANNEL (NADH, FAD), construct a side-by-side stack.
    % Otherwise just do the one channel (below)
    if ~isempty(Dnad)&~isempty(Dfad)
        % If acquisition was stopped part way through, just go up to
        % the last frame that had both NADH and FAD.
        sdtL = min([length(Dnad) length(Dfad)]);
        test = cell(sdtL,1);
        
        for i = 1:sdtL
            % Save time by skipping frames that were already done
            if isempty(find(frs2do==i))
                continue;
            end
            %             try
            test{i} = i;
            nsdt = bh_readsetup([sdtpath Dpos(posnum).name '\' Dnad(i).name]);
            fsdt = bh_readsetup([sdtpath Dpos(posnum).name '\' Dfad(i).name]);
            ve = nsdt.Version;
            dashes = strfind(Dnad(i).name,'_');
            t = str2num(Dnad(i).name(dashes(1)+1:dashes(2)-1));
            z = str2num(Dnad(i).name(dashes(3)+1:end-4));
            disp(['t' num2str(t) ', Pos' num2str(PosNum) ', z' num2str(z)]);
            
            % Find correct nameind index and fill in the new stack
            % frame number. Define inds in a cell so code can run
            % in parallel
            inds{i} = (PosInd)&([nameinds{:,2}]==t)&([nameinds{:,5}]==z);
            % Check Column 7 for '-1' flag. If it's not there, write the
            % tiff
            if ~isempty(find(cell2mat(nameinds(inds{i},7))==-1))
                continue;
            end
            % Only look at NADH image. If it's a fifo image, the corresponding FAD
            % better be, too.
            try % If bh read errors, likely the sdt was bad. Flag frame
                nadch1 = bh_getdatablock_v095(nsdt,block_1);
                if FLIMageBool nadFLIMage = QuickFLIMage(nadch1,1,15); end
                nadch1 = squeeze(sum(nadch1,1));
                
                
                fadch1 = bh_getdatablock_v095(fsdt,block_1);
                if FLIMageBool fadFLIMage = QuickFLIMage(fadch1,1,15); end
                fadch1 = squeeze(sum(fadch1,1));
                
                if SHGBool
                    SHGim = bh_getdatablock_v095(fsdt,block_2);
                    if bit16Bool
                        SHGim = uint16(squeeze(sum(SHGim,1)));
                    else
                        SHGim = uint8(squeeze(sum(SHGim,1)));
                    end
                end
            catch
                nameinds(inds{i},[7 8 9]) = num2cell(-1);
                continue;
            end
            % If either image is empty, skip this frame
            if isempty(find(nadch1))|isempty(find(fadch1))
                nameinds(inds{i},[7 8 9]) = num2cell(-1);
                continue;
            end
            
            if str2num(ve(1))==3
                if bit16Bool
                    im = uint16([(nadch1) (fadch1)]);
                else
                    im = uint8([(nadch1) (fadch1)]);
                end
                if FLIMageBool FLim = uint8([(nadFLIMage) (fadFLIMage)]); end
            else
                error('Data type not a fifo image');
            end
            imwrite(im,[tifpath 'fr' num2str(i,'%05i') '.tif'],'tiff','Compression','none');
            if FLIMageBool imwrite(FLim,[FLIMagepath 'fr' num2str(i,'%05i') '.tif'],'tiff','Compression','none'); end
            % Sometimes SHG images have literally, like 0 photons. If so,
            % set SHGim to a matrix of zeros to avoid an error.
            if SHGBool 
                if isempty(SHGim) SHGim = zeros(size(im,1),size(im,1)); end
                imwrite(SHGim,[SHGpath 'fr' num2str(i,'%05i') '.tif'],'tiff','Compression','none'); 
            end
            
            fclose('all');
            
        end
        
        % When tiff conversions are finally finished, fill in all the
        % nameinds and save below
        for i = 1:length(inds)
            if ~isempty(inds{i})
                indall = find(inds{i}); %Include flagged entries (-1)
                indnonflag = find(inds{i}&[nameinds{:,7}]>-1); %Exclude -1's
                Nindfin = indall(strcmp(nameinds(indall,4),'NADH'));
                Findfin = indall(strcmp(nameinds(indall,4),'FAD'));
                nameinds(indall,6) = num2cell(i);
                if size(indnonflag,2)==size(indall,2) % Neither NADH or FAD flagged
                    % Do nothing unless there's a flag
                else % One is flagged, so this frame is bad. Flag the other, too
                    nameinds(Nindfin,[7 8 9]) = num2cell(-1);
                    nameinds(Findfin,[7 8 9]) = num2cell(-1);
                end
            end
        end
    else % Assume 1 channel
        sdtL = max([length(Dnad) length(Dfad) length(Duser) length(D)]);
        for i = 1:sdtL
            if isempty(find(frs2do==i))
                continue;
            end
            if ~isempty(Dnad) % Only NADH images taken
                sdt = bh_readsetup([sdtpath Dpos(posnum).name '\' Dnad(i).name]);
                dashes = strfind(Dnad(i).name,'_');
                t = str2num(Dnad(i).name(dashes(1)+1:dashes(2)-1));
                z = str2num(Dnad(i).name(dashes(3)+1:end-4));
            end
            if ~isempty(Dfad)% Only FAD images taken
                sdt = bh_readsetup([sdtpath Dpos(posnum).name '\' Dfad(i).name]);
                dashes = strfind(Dfad(i).name,'_');
                t = str2num(Dfad(i).name(dashes(1)+1:dashes(2)-1));
                z = str2num(Dfad(i).name(dashes(3)+1:end-4));
            end
            if ~isempty(Duser)% Only FAD images taken
                sdt = bh_readsetup([sdtpath Dpos(posnum).name '\' Duser(i).name]);
                dashes = strfind(Duser(i).name,'_');
                t = str2num(Duser(i).name(dashes(1)+1:dashes(2)-1));
                z = str2num(Duser(i).name(dashes(3)+1:end-4));
            end
            ve = sdt.Version;
            
            disp(['t' num2str(t) ', Pos' PosNum ', z' num2str(z)]);
            % Find correct nameind index and fill in the new stack
            % frame number
            inds{i} = (PosInd)&([nameinds{:,2}]==t)&([nameinds{:,5}]==z);
            
            % Only look at NADH image. If it's a fifo image, the corresponding FAD
            % better be, too.
            ch1 = bh_getdatablock_v095(sdt,block_1);
            if FLIMageBool FLIMage = QuickFLIMage(ch1,1,15); end
            ch1 = squeeze(sum(ch1,1));
            % If either image is empty, skip this frame
            if isempty(find(ch1))
                nameinds(inds{i},[7 8 9]) = num2cell(-1);
                continue;
            end
            if SHGBool
                SHGim = bh_getdatablock_v095(sdt,block_2);
                if bit16Bool
                    SHGim = uint16(squeeze(sum(SHGim,1)));
                else
                    SHGim = uint8(squeeze(sum(SHGim,1)));
                end
            end
            if str2num(ve(1))==3
                if bit16Bool
                    im = uint16(ch1);
                else
                    im = uint8(ch1);
                end
                if FLIMageBool FLim = uint8(FLIMage); end
            else
                error('Data type not a fifo image');
            end
            imwrite(im,[tifpath 'fr' num2str(i,'%05i') '.tif'],'tiff','Compression','none');
            if FLIMageBool imwrite(FLim,[FLIMagepath 'fr' num2str(i,'%05i') '.tif'],'tiff','Compression','none'); end
            if SHGBool imwrite(SHGim,[SHGpath 'fr' num2str(i,'%05i') '.tif'],'tiff','Compression','none'); end;
            
            fclose('all');
        end
        % When tiff conversions are finally finished, fill in all the
        % nameinds and save below
        for i = 1:length(inds)
            if ~isempty(inds{i})
                indall = find(inds{i}); %Include flagged entries (-1)
                %                 indnonflag = find(inds{i}&[nameinds{:,7}]>-1); %Exclude -1's
                %                 Uindfin = indall(strcmp(nameinds(indall,4),'UserChan'));
                nameinds(indall,6) = num2cell(i);
                %                 if size(indnonflag,2)~=size(indall,2)
                %                 % If it was flagged...
                %                     nameinds(Uindfin,[7 8 9]) = num2cell(-1);
                %                 end
            end
        end
    end
end
% Any file not converted to a tiff will still have '-1' in column 6.
% Propagate that flage to columns 7 and 8 to exclude from future processing
nameinds([nameinds{:,6}]==-1,[7 8 9]) = num2cell(-1);
save([acqpath 'multiD_indices.mat'],'nameinds')

%% Z-binning stack
% Optional, make additional stacks in 'sorted_sdt' folder where z frames
% are binned together by adding, to produce one frame per time point
for posnum = 1:length(Dpos)
    % Get xdim and ydim
    clear StIms StIm tifpath Dtiffs ImInf xdim ydim PosNum PosInd Tind Zs
    PosNum = Dpos(posnum).name(4:end);
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    tifpath = [sdtpath 'IntTiffs_' Dpos(posnum).name '_' Run '\'];
    Dtiffs = dir(tifpath); Dtiffs(1:2)=[];
    SHGpath = [sdtpath 'SHGTiffs_' Dpos(posnum).name '_' Run '\'];
    
    ImInf = imfinfo([tifpath Dtiffs(1).name]);
    xdim = ImInf.Width; ydim = ImInf.Height;
    Tind = unique([nameinds{PosInd & [nameinds{:,7}]>-1, 2}]);
    WrtName = [UpOneDir([sdtpath 'IntTiffs_' Dpos(posnum).name '_' Run '\']) 'Pos' PosNum '_zBin.tif'];
    WrtNameSHG = [UpOneDir([sdtpath 'IntTiffs_' Dpos(posnum).name '_' Run '\']) 'SHG_Pos' PosNum '_zBin.tif'];
    
    for t = Tind
        Zs = unique([nameinds{PosInd & [nameinds{:,7}]>-1,5}]);
        StIm = zeros(ydim,xdim);
        StImSHG = zeros(ydim,ydim);
        for z = Zs
            FrNum = unique([nameinds{find(PosInd & [nameinds{:,2}]==t & [nameinds{:,5}]==z & [nameinds{:,7}]>-1),6}]);
            if ~isempty(FrNum)
                im = imread([tifpath 'fr' num2str(FrNum,'%05i') '.tif']);
                StIm = StIm + double(im);
                if SHGBool
                    im = imread([SHGpath 'fr' num2str(FrNum,'%05i') '.tif']);
                    StImSHG = StImSHG + double(im);
                end
            end
        end
        if bit16Bool
            StIms(t+1).data = uint16(StIm);
            if SHGBool StImsSHG(t+1).data = uint16(StImSHG); end
        else
            StIms(t+1).data = uint8(StIm);
            if SHGBool StImsSHG(t+1).data = uint8(StImSHG); end
        end
        
    end
    StkWrite(StIms, WrtName)
    if SHGBool StkWrite(StImsSHG, WrtNameSHG); end
end

%% Ruth dye
% Are there Ruthenium dye frames? If so, convert them independently
% (typically one middle z-plane take for each time point)
Druth = dir([sdtpath Dpos(posnum).name '\*Ruth*.sdt']);
if ~isempty(Druth)
    for posnum = 1:length(Dpos)
        Ruthpath = [sdtpath 'RuthTiffs_' Dpos(posnum).name '_' Run '\'];
        [a,b] = mkdir(Ruthpath);
        DRuthTiffs = dir([Ruthpath '\*.tif']);
        for i = 1:length(Druth)
            D = dir([Ruthpath 'T' num2str(i,'%05i') '.tif']);
            if isempty(D)
                Rsdt = bh_readsetup([sdtpath Dpos(posnum).name '\' Druth(i).name]);
                Ruthch1 = bh_getdatablock_v095(Rsdt,block_1);
                if bit16Bool
                    Ruthch1 = uint16(squeeze(sum(Ruthch1,1)));
                else
                    Ruthch1 = uint8(squeeze(sum(Ruthch1,1)));
                end
                imwrite(Ruthch1,[Ruthpath 'T' num2str(i,'%05i') '.tif'],'tiff','Compression','none');
            end
        end
    end
end

%% WoW Circular Crops
% Found it is useful for all analysis to go ahead and get the crops using
% image cross-correlation right up front. ImageJ Weka segmentation will use
% this info.
if WellCropBool
    NADHBool = 1; FADBool = 1;
    if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
    if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
    for posnum = 1:length(Dpos)
        % Find peripheral well area and exclude
        tifpath = [sdtpath 'IntTiffs_' Dpos(posnum).name '_' Run '\'];
        Dtiffs = dir(tifpath); Dtiffs(1:2)=[];
        
        im = double(imread([tifpath Dtiffs(1).name]));
        xdim1ch = size(im,2); ydim = size(im,1);
        if NADHBool&FADBool
            xdim1ch = xdim1ch/2;
            im1ch = im(:,xdim1ch+1:end);
        elseif xdim1ch==ydim
            xdim1ch = xdim1ch;
            im1ch = im;
        else
            error('Something wrong with image dimensions.')
        end
        [Wmask,acont,CoM,rad] = WoWcropActiveCont(im1ch);
        Cx = CoM(1); Cy = CoM(2);
        imwrite(Wmask,[sdtpath 'Crop_' Dpos(posnum).name '_' Run '.tif'],'tiff','Compression','none');
        save([acqpath 'WellCenters.mat'],'Cx','Cy');
    end
end




