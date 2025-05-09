function Ilastik_TrainImGen
% Just takes a specified image, selects the NADH or FAD channel, scales it
% by the IllProfCal in the data cal folder, and resaves it in the sample
% folder for use in training Ilastik

% % % Test settings:
% clear all;
% IntTiffsPath = 'C:\Dropbox\data\s1_a1_zscan\sorted_sdts\IntTiffs_Pos0_s1_a1_zscan\';
% fn = 'fr00015.tif';
% Nbool=0; Fbool=1; Pbool=0; Obool=0; NrmBool=1; SepMskBool=0;
% thresh = [100 100 100 100];

% User picks:
[fn,IntTiffsPath]=uigetfile('*.tif','Select the tiff(s) to create training images','MultiSelect','on');
[Nbool,Fbool,Pbool,Obool,NrmBool,SepMskBool,ThreshLow,ThreshHigh] = IllastikTrainingGUI;
thresh = [ThreshLow ThreshHigh ThreshLow ThreshHigh]; % same for both

if Nbool
    chan = 'NADH'; chnum = 1;
elseif Fbool
    chan = 'FAD'; chnum = 2;
elseif Pbool
    chan = 'Product'; chnum = [1 2];
elseif Obool
    chan = 'User'; chnum = 3;
end
if NrmBool
    NrmLab = '_nrm';
else
    NrmLab = '';
end
% chan = questdlg('NADH or FAD training','???NADH or FAD???','NADH','FAD','Product','NADH'); % Previous dialog box
if ~ischar(IntTiffsPath)|isempty(chan)
    disp('Images or channels not selected.');
    return;
end
Gblur=12;
G = fspecial('gaussian',[Gblur Gblur],Gblur);

if ~iscell(fn) fn={fn}; end

% Data paths
if IntTiffsPath(end)~='\'; IntTiffsPath = [IntTiffsPath '\']; end;
acqpath = UpOneDir(UpOneDir(IntTiffsPath));
slashes = strfind(IntTiffsPath,'\');
Lab = IntTiffsPath(slashes(end-1)+1:slashes(end)-1);
% Determine what kind of acquisition this is. MultiD, NADH, FAD, or just
% generic series
NADHBool = 1; FADBool = 1;
try
    load([acqpath 'multiD_indices.mat']);
catch
    nameinds = nan(1,9);
end
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if ~exist('thresh')|thresh==-1 thresh = [0 10^6 0 10^6]; end % Unless specified, brightest 2 clusters is pretty good at getting mitochondria and excluding nucleus
sdtpath = [acqpath 'sorted_sdts\'];

% Find cal_files folder
if exist([acqpath 'cal_files'])==7 % if cal_files is in acqpath
    calpath = [acqpath 'cal_files\'];
elseif exist([IntTiffsPath 'cal_files'])==7 % if cal_files in tiffs folder
    calpath = [IntTiffsPath 'cal_files\'];
elseif exist([UpOneDir(acqpath) 'cal_files'])==7 % if cal_files in daypath
    calpath = [UpOneDir(acqpath) 'cal_files\'];
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

% Automatically search for Crop tiffs that would have been created in
% MultiD_tiff_convert
PosInd = strfind(IntTiffsPath,'Pos');
dashes = strfind(IntTiffsPath,'_');
try % If not a MultiD, then assume no crop
nextdash = dashes(dashes>PosInd);
PosNum = IntTiffsPath(PosInd+3:nextdash(1)-1);
im0 = double(imread([IntTiffsPath fn{1}]));
ydim = size(im0,1);
D = dir([sdtpath '*Crop*Pos' PosNum '*']);
if ~isempty(D)
    Crop = imread([sdtpath D(1).name]);
else
    Crop = 255.*ones(ydim,ydim);
end
catch
    finf = imfinfo([IntTiffsPath '\' fn{1}]);
    ydim = finf.Height;
    Crop = 255.*ones(ydim,ydim); 
end

% Check for custom channel
chans = unique(nameinds(:,4));
chans2 = chans; chans2(strcmp(chans2,'NADH'))=[]; chans2(strcmp(chans2,'FAD'))=[];
UserChBool = 0; if ~isempty(chans2) UserChBool = 1; chnum =3; end

% Load the frames and scale them by IllProfCal, then save
for i = 1:length(fn)
    
    im0 = double(imread([IntTiffsPath fn{i}]));
    
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
    elseif Obool
        im(:,:,3) = im0;
    else
        error('Channel error')
    end
    
    % Look for number of scans to create images. If scan=-1, frame is a
    % dud, so just divide dud image by the av number of scans for the
    % acq.
    fr = str2num(fn{i}(3:end-4));
    if NADHBool nscans(1) = double(nameinds{[nameinds{:,6}]'==fr & strcmp(nameinds(:,3),PosNum) & strcmp(nameinds(:,4),'NADH'),9}); end
    if FADBool nscans(2) = double(nameinds{[nameinds{:,6}]'==fr & strcmp(nameinds(:,3),PosNum) & strcmp(nameinds(:,4),'FAD'),9}); end
    try
        if UserChBool nscans(3) = double(nameinds{[nameinds{:,6}]'==fr & strcmp(nameinds(:,3),PosNum) & strcmp(nameinds(:,4),chans2{1}),9}); end
    catch
        nscans = [1 1 1]; % If intensity series, make generic
        chan = 'User';
    end
    % Process all channels in a loop (easier to code), but below, just save
    % the one specified by the user.
    for ch = find([NADHBool FADBool UserChBool])
        
        % Scale, divide by number of scans.
        im(:,:,ch) = im(:,:,ch)./nscans(ch);
        
        % Crop well
        im(find(~Crop)+numel(im(:,:,1))*(ch-1))=0;
        
        % Only do selected channel
        if isempty(find(chnum==ch)) continue; end
        
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
            
            % Tried at some point, but probably not useful. Take out soon.
%             % Divide background by image mean to keep that relative to
%             % foreground signal
%             BGmsk = ~singlemask; BGind = find(BGmsk)+numel(im(:,:,1))*(ch-1);
%             im(BGind) = im(BGind)./mean(Mns);
        end
    end
    
    % Scale by IllProfCal file. Note, we're doing this after the
    % normalization and 1st-order segmentation because I found IllProfCal
    % division actually introduced artifacts for the bpfilter.
    if ProfBool
        for k = 1:size(im,3) im(:,:,k) = im(:,:,k)./IllProfCal.*mean(IllProfCal(:)); end
    end
    
    % Pick out unique day number at add it to training image file name.
    % This avoids overwriting previous training images with same file name.
    daypath = UpOneDir(UpOneDir(calpath));
    slashes = strfind(daypath,'\');
    dayname = daypath(slashes(end-1)+1:end-1);
    dashes = strfind(dayname,'_');
    dayind = strfind(dayname,'_D');
    daynum = dayname(dayind+2:dashes(end)-1);
    if strcmp(chan,'NADH')
        %         imwrite2tif(im(:,:,1),[],[acqpath '\IlTrain_' fn{i}(1:end-4) '_N_' Lab NrmLab '.tif'],'single')
        imwrite2tif(im(:,:,1),[],[calpath '\IlTrain_' fn{i}(1:end-4) '_D' daynum '_P' PosNum '_N' NrmLab '.tif'],'single')
    elseif strcmp(chan,'FAD')
        %         imwrite2tif(im(:,:,2),[],[acqpath '\IlTrain_' fn{i}(1:end-4) '_F_' Lab NrmLab '.tif'],'single')
        imwrite2tif(im(:,:,2),[],[calpath '\IlTrain_' fn{i}(1:end-4) '_D' daynum '_P' PosNum '_F' NrmLab '.tif'],'single')
    elseif strcmp(chan,'Product')
        %         imwrite2tif(im(:,:,1).*im(:,:,2),[],[acqpath '\IlTrain_' fn{i}(1:end-4) '_P_' Lab NrmLab '.tif'],'single')
        imwrite2tif(im(:,:,1).*im(:,:,2),[],[calpath '\IlTrain_' fn{i}(1:end-4) '_D' daynum '_P' PosNum '_P' NrmLab '.tif'],'single')
    elseif strcmp(chan,'User')
        %         imwrite2tif(im(:,:,1).*im(:,:,2),[],[acqpath '\IlTrain_' fn{i}(1:end-4) '_P_' Lab NrmLab '.tif'],'single')
        imwrite2tif(im(:,:,3),[],[calpath '\IlTrain_' fn{i}(1:end-4) '_D' daynum '_P' PosNum '_U' NrmLab '.tif'],'single')
    else
        error('Problem with channels... NOOOOO!!!');
    end
    
end

end % Function end


function SettingsCheckBox
% Create figure
h.f = figure('units','pixels','position',[200,200,150,50],...
    'toolbar','none','menu','none');
% Create yes/no checkboxes
h.c(1) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,30,50,15],'string','yes');
h.c(2) = uicontrol('style','checkbox','units','pixels',...
    'position',[90,30,50,15],'string','no');
% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[40,5,70,20],'string','OK',...
    'callback',@p_call);
% Pushbutton callback
    function p_call(varargin)
        vals = get(h.c,'Value');
        checked = find([vals{:}]);
        if isempty(checked)
            checked = 'none';
        end
        disp(checked)
    end
end



