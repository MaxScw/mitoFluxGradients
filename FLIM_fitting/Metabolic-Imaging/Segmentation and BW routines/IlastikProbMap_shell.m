function IlastikProbMap_shell(acqpath,ModName,poss,Frs,thresh,SepMskBool)
% Ilastik (http://ilastik.org) is a freeware analysis suite for pixel segmentation.
% This routine comes after tif creation from sdt's. It calls Ilastik via
% the shell from within Matlab. Ilastik calculates probability maps that
% pixels are mitochondria, cytoplasm, or background. The maps are saved to
% the data folder.
% Inputs:
%  -acqpath: path to folder
%  -ModName: path for the ilastik 'project' file. This specifies the trained
%         model used to segment the images.
%         Note: if a model with ModName is in the cal_files directory, that
%         one is chosen. Otherwise, program looks in the default Common dir.
%  -Frs: Maybe you want to do only a range of frames
%  -thresh: upper and lower threshold values for NADH and FAD [Nmin Nmax Fmin Fmax]
%  -SepMskBool: If images contained multiple masks that whose intensities
%         were normalized separately, enter '1'

% % TEST in script mode:
% % clear all;
% acqpath = 'C:\Dropbox\data\s1_a1_zscan\';
% ModName = 'Mod_N_nrm';
% Frs = [15];
% poss= 0;

% Path to ilastik executable
if exist('C:\Program Files\ilastik-1.3.2post1\')==7
    IlCom = '"C:\Program Files\ilastik-1.3.2post1\ilastik.exe" --headless ';
else
    IlCom = '"C:\Program Files\ilastik-1.3.0\ilastik.exe" --headless ';
end
ModPath = 'Z:\Lab\Common\Metabolic Imaging Calibration Files\Ilastik_files\';
% ModPath = 'C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\IlastikMods\';

% Data paths
if acqpath(end)~='\'; acqpath = [acqpath '\']; end;
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end
NADHBool = 1; FADBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if ~exist('area_cuts')|area_cuts==-1 area_cuts = [1000 10^6]; end
if ~exist('thresh')|thresh==-1 thresh = [100 100 100 100]; end % Units of std dev
if ~exist('SepMskBool')|SepMskBool==-1 SepMskBool = 0; end % Units of std dev
slashes = strfind(acqpath,'\');
Run = acqpath(slashes (end-1)+1:end-1);
sdtpath = [acqpath 'sorted_sdts\'];
Gblur = 12;
G = fspecial('gaussian',[Gblur Gblur],Gblur); % Blur function 

% Find cal_files folder
if exist([acqpath 'cal_files'])==7 % if cal_files is in acqpath
    calpath = [acqpath 'cal_files\'];
elseif exist([UpOneDir(acqpath) 'cal_files'])==7 % if cal_files in daypath
    calpath  = [UpOneDir(acqpath) 'cal_files\'];
else
    error('Cannot locate cal_files folder. Place in daypath or acqpath, please');
end

if ~exist('ModName')|ModName==-1 error('You really should specify a Model'); end;
if ~isempty(strfind(ModName,'nrm')) NrmBool=1; else NrmBool=0; end;
if isempty(strfind(ModName,'.ilp')) ModName = [ModName '.ilp']; end;
if exist([ModPath ModName],'file')==0 ModPath = 'C:\Ilastik_files\'; end % Work Comp vs home comp check
if exist([calpath ModName],'file')==2 ModPath = calpath; end % cal_files dir check

% Load IllProfCal.mat file
ProfBool = 0;
Dprof = dir([calpath 'IllProfCal*.mat']);
if ~isempty(Dprof)
    load([calpath Dprof(1).name]);
    IllProfCal = double(IllProfCal);
    %         IllProfCaldual = [IllProfCal IllProfCal];
    ProfBool = 1;
end

% Determine channel
if strfind(ModName,'_N')
    chan = 'NADH';
elseif strfind(ModName,'_F')
    chan = 'FAD';
elseif strfind(ModName,'_P')
    chan = 'Product';
else
    chan = 'User';
end
ModName = [' --project="' ModPath ModName '"']


Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    if ~strcmp(Dpos(i).name(1:3),'Pos')
        remove = [remove i];
    end
end
Dpos(remove)=[];

if ~isempty(Dpos)
    srtpath = [acqpath 'sorted_sdts\'];
    for i = 1:size(Dpos,1)
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
        Probpaths{i} = [srtpath 'ProbMaps_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir([Probpaths{i}]);
        
        %         % save figures with overlaid ROIs to do quick checks after batch processing
        %         ROIpaths{i} = [srtpath 'ROIsCheckGen_' Dpos(i).name '_' Run '\'];
        %         [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = path;
end


% Check for custom channel
chans = unique(nameinds(:,4));
chans2 = chans; chans2(strcmp(chans2,'NADH'))=[]; chans2(strcmp(chans2,'FAD'))=[];
UserChBool = 0; if ~isempty(chans2) UserChBool = 1; end

for posnum = 1:size(Dpos,1)
    %     posnum =1;
    PosNum = Dpos(posnum).name(4:end);
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        for ps = 1:length(poss) PsBools(ps) = strcmp(num2str(poss(ps)),PosNum); end
        if isempty(find(PsBools))
            continue;
        end
    end
    
    frames = unique(sort([nameinds{PosInd&[nameinds{:,7}]>-1,6}]));
    frames = frames(frames>0);
    ts = unique(sort([nameinds{PosInd,2}]))+1;
    Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if exist('Frs')&Frs~=-1 frames = Frs; end;
    
    set(gcf,'PaperPositionMode', 'auto')
    clear Nmasks Fmasks Masks EggVals NADim FADim
    
    cd(acqpath);
    im1=imread([tifpaths{posnum} 'fr' num2str(frames(1),'%05i') '.tif']);
    [ydim,xdim] = size(im1);
    % Check for well crop file for this Pos
    D = dir([sdtpath '*Crop*Pos' PosNum '*']);
    if ~isempty(D)
        Crop = imread([sdtpath D(1).name]);
    else
        Crop = 255.*ones(ydim,ydim);
    end
    meannumscs = mean([nameinds{[nameinds{:,9}]>0,9}]);
        
    %% Find the frames, make temp copies of them, put their file names in a long string for ilastik call
    IlFiles = []; DelFiles = {};
    for i = 1:length(frames)
        t(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),2}]))+1;
        Z(i) = unique(sort([nameinds{PosInd&[nameinds{:,6}]==frames(i),5}]))+1;
        
        disp(['Pos' PosNum ', Fr ' num2str(frames(i))])
        
        im0 = double(imread([tifpaths{posnum} 'fr' num2str(frames(i),'%05i') '.tif']));
        % Reshape tiff image to load channels into the 3rd dim with
        % consistent indexing - NADH=1, FAD=2, CustomCh=3
        % Also get number of scans
        xdim = size(im0,2); ydim = size(im0,1);
        if NADHBool&FADBool
            im(:,:,1) = im0(:,1:ydim);
            im(:,:,2) = im0(:,ydim+1:2*ydim);
        elseif NADHBool & ~FADBool
            im(:,:,1) = im0;
        elseif ~NADHBool & FADBool
            im(:,:,2) = im0;
        elseif ~NADHBool & ~FADBool
            im(:,:,3) = im0;
        end
        
        % Look for number of scans to create images. If scan=-1, frame is a
        % dud, so just divide dud image by the av number of scans for the
        % acq.
        if NADHBool nscans(1) = double(nameinds{[nameinds{:,6}]'==frames(i) & strcmp(nameinds(:,3),PosNum) & strcmp(nameinds(:,4),'NADH'),9}); end
        if FADBool nscans(2) = double(nameinds{[nameinds{:,6}]'==frames(i) & strcmp(nameinds(:,3),PosNum) & strcmp(nameinds(:,4),'FAD'),9}); end
        if UserChBool nscans(3) = double(nameinds{[nameinds{:,6}]'==frames(i) & strcmp(nameinds(:,3),PosNum) & strcmp(nameinds(:,4),chans2{1}),9}); end
        % Process all channels in a loop (easier to code), but below, just save
        % the one specified by the user.
        for ch = find([NADHBool FADBool UserChBool])
            % Scale, divide by number of scans.
            im(:,:,ch) = im(:,:,ch)./nscans(ch);
            
            % Crop well
            im(find(~Crop)+numel(im(:,:,1))*(ch-1))=0;
            
            % Get thresh pix. Run simple gauss blur first so that we're
            % removing features, not fluctuations (e.g. polar body)
            gim = imfilter(im(:,:,ch),G,'same');
            close all;
            
            % Normalize fluorescence intensity distribution?
            if NrmBool
                % Do a basic segmentation separate embryos from background
                im2=im(:,:,ch); im2 = im2./mean(im2(:)); % Initial scale to mean
                [masks singlemask num] = BpassBGSeg(im2,-1,-1,0);
                % Note the 'ch-1' is because NADH images tend to need a
                % lower threshold because we want to include more
                % cytoplasm.
%             % Previously used kmeans, but kind of unstable
%             [masks singlemask num] = Masks_Kmeans_FLIMages(im(:,:,ch),-1,Gblur,3,[2 3],[100 10^8]);
%             % Quick check, if kmeans got more than 70% of the FOV, it was
%             % probably too agressive, so repeat with 2 groups instead of 3
%             if length(find(singlemask))>length(find(gim))*.7
%                 [masks singlemask num] = Masks_Kmeans_FLIMages(im(:,:,ch),-1,Gblur,2,2,[100 10^8]);
%             end
                
                % Scale to mean of fluorpix. (find(bw))
                % Normalize each individual mask's int dist to 1
                PrevDilMsks = zeros(size(im,1),size(im,2)); % Joint mask for accounting for overlap of dilated masks
                
                if ~SepMskBool clear masks; masks{1} = singlemask; end
                
                for m = 1:length(masks)
                    msk = imdilate(masks{m},strel('disk',8)); % Dilate the individual mask
                    %                 subplot(1,3,1); imshow(msk,[]);
                    msk(find(msk&PrevDilMsks))=0;
                    %                 subplot(1,3,2); imshow(msk,[]);
                    PrevDilMsks = PrevDilMsks | msk;
                    %                 subplot(1,3,3); imshow(PrevDilMsks,[]);
                    
                    % Apply thresholds (units of std dev)
                    mskmn = mean(gim(find(msk))); mskstd = std(gim(find(msk)));
                    %                 subplot(1,2,1); imshow(im(:,:,ch),[]);
                    threshpix = intersect( find(gim<mskmn-mskstd*thresh(1) | (gim>mskmn+mskstd*thresh(2))) , find(msk) );
                    im(threshpix+numel(im(:,:,1))*(ch-1))=0;
                    %                 subplot(1,2,2); imshow(im(:,:,ch),[]);
                    
                    % Normalize to mean
                    mind = find(msk)+numel(im(:,:,1))*(ch-1);
                    %                 Mns(m) = mean(im(mind));
                    im(mind) = im(mind)./mean(im(mind));
                end
                % Clear everything that fell into kmeans background. BG
                % pixels can sometimes get confused for dim cyto signal
                bind = find(~singlemask)+numel(im(:,:,1))*(ch-1);
                im(bind) = 0;
            end
        end
        
        % Scale by IllProfCal file. Note, we're doing this after the
        % normalization and 1st-order segmentation because I found IllProfCal
        % division actually introduced artifacts for the bpfilter.
        if ProfBool
            for k = 1:size(im,3) im(:,:,k) = im(:,:,k)./IllProfCal.*mean(IllProfCal(:)); end
        end
        
        cd(tifpaths{posnum});
        if strcmp(chan,'NADH')
            imwrite2tif(im(:,:,1),[],['fr' num2str(frames(i),'%05i') '_N.tif'],'single')
            IlFiles = [IlFiles ' fr' num2str(frames(i),'%05i') '_N.tif'];
            DelFiles = [DelFiles; 'fr' num2str(frames(i),'%05i') '_N.tif'];
        elseif strcmp(chan,'FAD')
            imwrite2tif(im(:,:,2),[],['fr' num2str(frames(i),'%05i') '_F.tif'],'single')
            IlFiles = [IlFiles ' fr' num2str(frames(i),'%05i') '_F.tif'];
            DelFiles = [DelFiles; 'fr' num2str(frames(i),'%05i') '_F.tif'];
        elseif strcmp(chan,'Product')
            imwrite2tif(im(:,:,1).*im(:,:,2),[],['fr' num2str(frames(i),'%05i') '_P.tif'],'single')
            IlFiles = [IlFiles ' fr' num2str(frames(i),'%05i') '_P.tif'];
            DelFiles = [DelFiles; 'fr' num2str(frames(i),'%05i') '_P.tif'];
        else % Other, user channel
            imwrite2tif(im(:,:,3),[],['fr' num2str(frames(i),'%05i') '_U.tif'],'single')
            IlFiles = [IlFiles ' fr' num2str(frames(i),'%05i') '_U.tif'];
            DelFiles = [DelFiles; 'fr' num2str(frames(i),'%05i') '_U.tif'];
        end
    end
    
    % Call Ilastik to make the prob maps
    status = dos([IlCom ModName IlFiles]);
    
    % Delete the temp files in the folder and move the prob maps
    for i = 1:size(DelFiles,1)
        delete(DelFiles{i});
        movefile([DelFiles{i}(1:end-4) '_Probabilities.tif'],[Probpaths{posnum} DelFiles{i}(1:end-4) '_Probabilities.tif']);
    end
end



