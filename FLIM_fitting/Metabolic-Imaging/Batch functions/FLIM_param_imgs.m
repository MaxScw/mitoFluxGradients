% function FLIM_param_imgs(acqpath,poss,MaskLab,Bin,PhotLim,fxshftBool,fxBgBool,nexpo,TplotBool)

% Incomplete, but getting there. I completed the pixel-by-pixel fitting,
% but it takes forever and it also gave weird FLIMage parameter results.

% INPUTs:
% -acqpath: path to acquisition
% -poss: path to acquisition
% -MaskLab: path to acquisition
% -MinFrames: path to acquisition
% -ProbThreshMeth - 1=prob map, 2=prob^2, any number between 0 and 1
%   specifies a custom prob threshold. DEFAULT value is 0.7 hard thresh.
% -MasksToDo - cell of arrays of mask numbers to get decays for. Each cell
%   element has an array of mask numbers for the corresponding position. For
%   a given position, if you want to do all the masks, enter '-1'.
%   NOTE: array must be in the same order as the '..._Masks.mat' files in the
%   acquisition acqpath. For example, if there are Pos0, Pos1, and Pos2, but
%   Pos0 didn't have any data, you would order a 2-element cell corresponding
%   to {[Pos1MasksToDo],[Pos2MasksToDo]}

% TEST:
clear all;
acqpath = 'C:\Dropbox\data\s1_a1';
MaskLab = 'JointMasks';
Bin=3;
PhotLim=15;
MaxIter=10;
% MasksToDo = {[-1],[2 4],[],[],[],[],[]};

try parpool; catch; end

if acqpath(end)~='\' acqpath = [acqpath '\']; end;
slashes = strfind(acqpath,'\');
Run = acqpath(slashes(end-1)+1:end-1);
if ~exist('MaskLab')|MaskLab==-1 MaskLab = 'Masks'; end
if ~exist('MinFrames')|MinFrames==-1 MinFrames = 1; end
if ~exist('ProbThreshMeth')|ProbThreshMeth==-1 ProbThreshMeth = .7; end
if ~exist('fxBools')|fxBools==-1 fxBools=[1; zeros(10,1)]; end % Default, fix only shift
IRfact=1;

% Load nameinds to be sure to connect correct mask to each sdt file
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end

% Determine fix parameters, load or get values
if ~isempty(find(fxBools==1))
    if ~exist([UpOneDir(acqpath) 'FixedParsVals.mat']) FLIMGetFixedPars(acqpath); end
    load([UpOneDir(acqpath) 'FixedParsVals.mat']);
end

% Find IRF sdt file for the day, convert to 'FLIMirf.mat'
% Also get bounds needed for fitting
% Find cal_files folder
if exist([acqpath 'cal_files'])==7 % if cal_files is in acqpath
    calpath = [acqpath 'cal_files\'];
elseif exist([UpOneDir(acqpath) 'cal_files'])==7 % if cal_files in daypath
    calpath  = [UpOneDir(acqpath) 'cal_files\'];
else
    error('Cannot locate cal_files folder. Place in daypath or acqpath, please');
end

% Determine which channels are present
uchans = unique(nameinds(:,4));
if ~isempty(find(strcmp(uchans,'NADH'))) Chind(1) = 1; end
if ~isempty(find(strcmp(uchans,'FAD'))) Chind(2) = 2; end
if ~isempty(find(strcmp(uchans,'UserChan'))) Chind(3) = 3; end
Chind(Chind==0)=[];
ChLabs = {'NADH','FAD','UserChan'};

% Load initial decay
D = dir([acqpath '*decays*.mat']);
load([acqpath D(1).name]);
decind = ~cellfun('isempty',decay_struct);decind = find(decind); decind = decind(1);
dectime = decay_struct{decind}.time;
dt = dectime(2)-dectime(1);
% Also get laser rep rate from any of the sdts in the folder
Dsdt = subdir([acqpath '\*.sdt']);
[bla,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup(Dsdt(1).name));
Dirfsdt = dir([calpath '*irf_*.sdt']);

% IRF
% Assume only 1 IRF in cal_files
if length(Dirfsdt)>1 warning('Should have only 1 IRF per acquisition! Wassamaddawichu?!'); end
if ~isempty(Dirfsdt)
    IRFstruct = FLIMLoadIRF([calpath Dirfsdt(1).name]);
    irftime = IRFstruct.time;
    [dt_irf,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup([calpath Dirfsdt(1).name])); % Assume the same even if using different IRFs
    % Calc adcratio
    adcratio = length(irftime)/length(dectime);
else
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

sdtpath = [acqpath 'sorted_sdts\'];
Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    % Check that folder is a 'Pos', and see if there is a corresponding 'Mask' file
    if ~strcmp(Dpos(i).name(1:3),'Pos')|isempty(dir([acqpath MaskLab '_' Dpos(i).name  '.mat']))
        remove = [remove i];
    end
end
Dpos(remove)=[];

% FIT INPUT PARAMETERS
% param(1) : shift of decay model from IRF (usually ranges from -10 to 10)
% Functions
% 1-expo: P = A*exp(-t/tau1) + (1-A)
% P = param(2)*exp(-t/param(3))+(1-param(2))
% param(2): fractional amp of expo decay(s). Fraction of photons coming
%           from fluorescence signal vs dark noise background. (1-A) = bg
% param(3): lifetime of the first decay
%
% 2-expo: P = A*[(1-f2)*exp(-t/tau1)+ f2*exp(-t/tau2)] + (1-A)
% P = param(2)*( (1-param(4))*exp(-t/param(3))+param(4)*exp(-t/(param(5))) )+(1-param(2));
% param(4) : fraction of the second decay
% param(5) : lifetime of second decay
%
% 3-expo: P = A(1-f2-f3)*exp(-t/tau1)+ A*f2*exp(-t/tau2) + A*f3*exp(-t/tau3) + (1-A)
% P = param(2)*( (1-param(4)-param(6))*exp(-t/param(3))+param(4)*exp(-t/(param(5)))+param(6)*exp(-t/(param(7))) )+(1-param(2));
% param(6) : fraction of the third decay
% param(7) : lifetime of third decay
% NOTE: exponentials are sorted in order of increasing lifetimes
nexpo = 2; % Change later if you want
N_param = 2*nexpo+1;
% If no IRF present, assume synthetic IRF and add a fit param for gaussian width
if isempty(Dirfsdt) N_param = N_param + 1; end
if nexpo == 1
    %         [shf,     A,    tau1]
    p_init0 = [0,       1,       1]';
    p_min0 =  [-10,   0.7,    0.05]';
    p_max0 =  [10,      1,       5]';
elseif nexpo == 2
    %         [shf,     A,    tau1,      f2,    tau2]
    p_init0 = [  0,     1,      .2,     0.3,       3]';
    p_min0 =  [-400,  0.7,    0.01,    0.01,    0.05]';
    p_max0 =  [400,     1,       1,       1,      10]';
elseif nexpo == 3
    %         [ shf,    A,    tau1,      f2,    tau2,     f3,    tau3]
    p_init0 = [   0,    1,      .2,     0.3,     0.5,     .3,       3]';
    p_min0 =  [-200,  0.7,    0.05,    0.01,       1,   0.01,       1]';
    p_max0 =  [ 200,    1,      .1,       1,       5,      1,      10]';
else
    errordlg('not available yet')
end

% Synthetic IRF?
if exist('IRFstruct')
    irf = IRFstruct.irf;
    IRFdt_ns = dt_irf;
    laserT_ns = LaserT;
    DecMax = -1;
else
    irf = -1; dt_irf = -1;
    % For synthetic IRF, fill in the last fit parameter as gaussian width
    % Initial value for IRF width measured empirically from a
    % typical IRF. ~90 bins on a 4096 bin array
    p_init0(N_param) = 100;
    p_min0(N_param) = 1;
    p_max0(N_param) = 300;
    
    % If no IRF present, lm_decay creates a synthetic gaussian IRF,
    % but it helps to do an initial centering by finding the knee
    sumdec  = sum(decay_struct{decind}.decay,2);
    tmp = find(sumdec==max(sumdec)); % offset empirical
    DecMax = tmp(1);
end

for posnum = 1:size(Dpos,1)
    
    clear decay_structs
    uManPos = Dpos(posnum).name(4:end);
    %     strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
    %     uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    PosInd = strcmp(nameinds(:,3),uManPos); PosInd = PosInd';
    
    % Make a separate folder for each param, also for each param std
    [a,b] = mkdir([sdtpath 'FLIMages_Pos' uManPos '\frc2']);
    [a,b] = mkdir([sdtpath 'FLIMages_Pos' uManPos '\LT1']);
    [a,b] = mkdir([sdtpath 'FLIMages_Pos' uManPos '\LT2']);
    [a,b] = mkdir([sdtpath 'FLIMages_Pos' uManPos '\std_frc2']);
    [a,b] = mkdir([sdtpath 'FLIMages_Pos' uManPos '\std_LT1']);
    [a,b] = mkdir([sdtpath 'FLIMages_Pos' uManPos '\std_LT2']);
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),uManPos)
            continue;
        end
    end
    
    % Load Masks - only use 'JointMasks' for this because we're analyzing
    % the whole image. We don't need individual masks for each embryo.
    load([acqpath 'JointMasks_' Dpos(posnum).name '.mat']);
    
    
    %% Load all data
    Dsdt = dir([acqpath 'sorted_sdts\' Dpos(posnum).name '\*.sdt']);
    NumOfStds = length(Dsdt);
    filenames = [];
    clear dashes t chans z StkFr
    for i = 1:NumOfStds
        filenames{i}=Dsdt(i).name;
        dashes = strfind(filenames{i},'_');
        t(i) = str2num(filenames{i}(dashes(1)+1:dashes(2)-1));
        chans{i} = filenames{i}(dashes(2)+1:dashes(3)-1);
        % Will be the same indexing for all msks
        z(i) = str2num(filenames{i}(dashes(3)+1:dashes(3)+4));
        
        % Use these and nameinds to get the stack slice number, which
        % corresponds with the 'Masks' structure elements, used below
        subset = find((PosInd)&([nameinds{:,2}]==t(i))&([nameinds{:,5}]==z(i))&([nameinds{:,7}]>-1));
        % Should be only 2 frames. If 1 or less, a frame got dropped, so
        % just set decay_structs{i,msk} to [];
        if isempty(strfind([nameinds{subset,4}],chans{i}))
            StkFr(i) = -1; % flag
        elseif strcmp(chans{i},'Ruth')
            StkFr(i) = nameinds{i,6};
            timestp(i) = nameinds{i,8};
        else
            if strfind(nameinds{subset(1),4},chans{i})
                StkFr(i) = nameinds{subset(1),6};
                timestp(i) = nameinds{subset(1),8};
            elseif strfind(nameinds{subset(2),4},chans{i})
                StkFr(i) = nameinds{subset(2),6};
                timestp(i) = nameinds{subset(2),8};
            else
                error('Something wrong with the name indexes')
            end
        end
    end
    uniqchans = unique(chans);
    
    %Image structure that contains all the information about the newly
    %opened images
    nummsks = 0;
    for i = 1:size(Masks,2) nummsks(i) = size(Masks{i},2); end
    nummsks = max(nummsks);
    decay_structs = cell(NumOfStds,nummsks);
    
    %Two photon image block
    block=1; %1:2pf, 2:SHG
    pathpar = [acqpath 'sorted_sdts\' Dpos(posnum).name '\'];
    filenames = filenames;
    for i = 1:NumOfStds
        % use the StkFr, -1 flag here to skip
        % frames that dropped a frame or had something go wrong.
        if StkFr(i)<0
            continue;
        end
        %if StkFr(i)==2
        %     disp('')
        %end
        disp([Run ', Pos' uManPos ', Frame ' num2str(StkFr(i))]);
        %load sdt file
        sdt = bh_readsetup([pathpar filenames{i}]);
        AcqRng = GetPhotonCollectTRange([pathpar filenames{i}]);
        
        % Get indexes
        dashes = strfind(filenames{i},'_');
        
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
        
        fit_region  = [round(resol*(1-meas.tac_lh/100)+1) round(resol*(1-meas.tac_ll/100)-1)];
        noise_region = [round(resol*(1-meas.tac_lh/100)+1) round(resol*(1-meas.tac_lh/100)+1)+5];
        
        %%
        %Update image structure
        image = img;
        filename = filenames{i};
        %         acqpath = acqpath;
        dt = dt;
        %image plot handle
        image_handle = -99;
        %1 if pixel selected, 0 if not
        selected_pixel = zeros(size(img));
        %Handles for plot showing selected pixels
        selected_pixel_handle = [];
        %FLIM Data
        flim = uint8(flim);
        
        %% Get decay for image from pixels contained in the mask.
        
        % Do for each msk. For single generic mask for whole pictures (eg
        % cumulus cells), store mask in first element of 'Masks' cell.
        if exist('MasksToDo')&(MasksToDo{posnum}~=-1)
            MaskInd = MasksToDo{posnum};
        else
            MaskInd = 1:size(Masks{StkFr(i)},2);
        end
        
        %         for msk = MaskInd
        msk=1;
        %time axis
        decay = zeros(size(flim,1),1);
        time = (1:length(decay))'*dt;
        
        % If this msk wasn't found for this frame, continue
        if msk>size(Masks{StkFr(i)},2)
            continue;
        end
        fn = fieldnames(Masks{StkFr(i)}(msk));
        if isempty(eval(['Masks{StkFr(i)}(msk).' fn{1}]))
            continue;
        end
        
        % Get mask
        selected_pixel=Masks{StkFr(i)}(msk).L;
        
        % Determine channel
        if strcmp(chans{i},'NADH')
            chnum = 1;
        elseif strcmp(chans{i},'FAD')
            chnum = 2;
        elseif strcmp(chans{i},'UserChan')
            chnum = 3;
        end
        
        
        % DECAY EXTRACTION
        %         imshow(squeeze(sum(flim,1)),[]);
        flimclr = permute(repmat(~selected_pixel,1,1,size(flim,1)),[3 1 2]);
        flim(flimclr)=0;
        %         figure;imshow(squeeze(sum(flim,1)),[]);
        [px1,px2] = find(selected_pixel);
        
        % BIN photons (blur) to improve decays and
        flimbin = zeros(size(flim));
        for p = 1:length(px1)
            % adjusted in case pixel is at the edge of the image
            px1adj = [max([1 (px1(p)-Bin)]) min([size(flim,2) (px1(p)+Bin)])];
            px2adj = [max([1 (px2(p)-Bin)]) min([size(flim,2) (px2(p)+Bin)])];
            flimbin(:,px1(p),px2(p)) = sum(sum(flim(:,px1adj(1):px1adj(2),px2adj(1):px2adj(2)),2),3);
        end
        %         % Test:
        %         close all; imshow(squeeze(sum(flim,1)),[]); figure; imshow(squeeze(sum(flimbin,1)),[]);
        
        
        fit_result = zeros(N_param,3);
        
        % Parameters to fix (shift is main one). Loop over fix options.
        for fx = 1:length(fxBools)
            if fxBools(fx)
                fit_result(fx,1) = fxvals(fx,chnum);
                fit_result(fx,3) = 1;
            end
        end
        
        free= ~fit_result(1:N_param,3);
        for i = find(free==0)
            p_init(i) = fit_result(i,1);
            p_max(i) = p_init(i)+1;
            p_min(i) = p_init(i)-1;
        end
        
        dp = zeros(N_param,1);
        dp(2:N_param) = 0.001*free(2:N_param);
        dp(1) = 1*free(1);
        if ~exist('IRFstruct') dp(end) = 1*free(end); end
        
        % FITS
        parfor p = 1:length(px1)
            decay = squeeze(flimbin(:,px1(p),px2(p)));
            if sum(decay)>PhotLim
                
                %weight on residual
                %weight = (fit_end-fit_start+1)/sqrt(decay(fit_start:fit_end)'*decay(fit_start:fit_end));
                nonzero_decay = decay;
                nonzero_decay(find(decay==0))=1;
                sigy = sqrt(nonzero_decay);
                weight = 1./sigy(fit_region(1):fit_region(2),:);
                
                % rough estimate of noise in noise region
                est_noise = mean(decay(noise_region(1):noise_region(2),1));
                
                % Update the noise initial estimate, use inputs above for the
                % rest of the parameters.
                p_init = p_init0; p_init(2) = 1-est_noise/max(decay(:,1));
                p_min = p_min0;
                p_max = p_max0;
                
                counts = sum(decay);
                
                [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst, converged] = ...
                    lm(@lm_decay_model,p_init,time,decay,weight,dp,p_min,p_max,[nexpo,counts,fit_region(1),fit_region(2),dt_irf,LaserT,DecMax],fit_region(1),fit_region(2),irf,MaxIter,0);
                
                p_fits{p} = p_fit;
                sigmas{p} = sigma_p;
                
                %                 if mod(round(length(px1)/100),p)==0
                %                     disp([num2str(round(p/length(px1)*100)) '% complete']);
                %                 end
                
                
                % Test NumIters with this code
                % for i = 1:20
                % [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst, converged] = lm(@lm_decay_model,p_init,time,decay,weight,dp,p_min,p_max,[nexpo,counts,fit_region(1),fit_region(2),dt_irf,LaserT,DecMax],fit_region(1),fit_region(2),irf,i*5,0);
                % pfits(:,:,i) = p_fit;
                % end
                % %
                % close all; plot(time,decay); hold on; plot(irftime,irf);
                % figure;
                % plot(squeeze(pfits(4,1,:)))
                
                
            end
        end
        
        % Fill in flimage
        flimage = zeros([size(flim,2) size(flim,3) 3]);
        sigmage = zeros([size(flim,2) size(flim,3) 3]);
        for p = 1:length(p_fits)
            if ~isempty(p_fits{p}) flimage(px1(p),px2(p),:) = p_fits{p}(3:5); end
            if ~isempty(sigmas{p}) sigmage(px1(p),px2(p),:) = sigmas{p}(3:5); end
        end
        imwrite2tif(flimage(:,:,1),[],[sdtpath 'FLIMages_Pos' uManPos '\LT1\' filename(1:end-3) 'tif'],'single');
        imwrite2tif(flimage(:,:,2),[],[sdtpath 'FLIMages_Pos' uManPos '\frc2\' filename(1:end-3) 'tif'],'single');
        imwrite2tif(flimage(:,:,3),[],[sdtpath 'FLIMages_Pos' uManPos '\LT2\' filename(1:end-3) 'tif'],'single');
        imwrite2tif(sigmage(:,:,1),[],[sdtpath 'FLIMages_Pos' uManPos '\std_LT1\' filename(1:end-3) 'tif'],'single');
        imwrite2tif(sigmage(:,:,2),[],[sdtpath 'FLIMages_Pos' uManPos '\std_frc2\' filename(1:end-3) 'tif'],'single');
        imwrite2tif(sigmage(:,:,3),[],[sdtpath 'FLIMages_Pos' uManPos '\std_LT2\' filename(1:end-3) 'tif'],'single');
        
        % Histograms
        
    end
end
