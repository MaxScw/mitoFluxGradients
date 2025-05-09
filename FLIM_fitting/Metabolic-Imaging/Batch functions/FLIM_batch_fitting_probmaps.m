function FLIM_batch_fitting_probmaps(acqpath,poss,MaskLab,Bin,fxBools,nexpo,TplotBool)
% Simple. Search acqpath for .mat files with 'decays' in them, perform fits
% for all contained decays and save as 'fits_...'
% Take guts from 'FittingGUI_ver3_4'
% INPUTS
% -poss: optional, fit for only one position.
% -MaskLab: Optional mask label. E.g. 'JointMasks' for one analyzing all
%   blobs in a given frame together as a single mask. Default is analyzing
%   each blob/mask as separate particles.
% -Bin: Bins together every 'Bin' frames into  a single decay and then
%   fits. Useful if photon count is low.
% -fxBools: array of bools for parameter fix options. Same order as
%   parameters, e.g. fix shift and long lifetime is [1 0 0 0 1];
% -nexpo: number of exponents to fit. 1, 2, and 3 are currently coded (in lm_decay.m)
% -TplotBool: Routine outputs time plots by default, but enter 0 if you don't want them.

% % TEST in script mode:
% clear all;
% acqpath = 'Z:\Lab\Tim\Boston IVF\Discarded_study\2019-07-05 40X zscans, 20X overnight\s1_a1_zscan';
% % Bin = 1;
% % fxBools = [1 0 0 0 0 1];
% MaskLab='JointMasks';

if acqpath(end)~='\' acqpath = [acqpath '\']; end
if ~exist('fxBools')|fxBools==-1 fxBools=[1; zeros(10,1)]; end % Default, fix only shift
if ~exist('Bin')|Bin==-1 Bin=1; end
if ~exist('MaskLab')|MaskLab==-1 MaskLab = 'JointMasks'; end 
if ~exist('MaskLab')|MaskLab==-1 MaskLab='mask'; end 
if ~exist('nexpo')|nexpo==-1 nexpo=2; end
if ~exist('TplotBool')|TplotBool==-1 TplotBool=1; end
% nexpo = 2; % Change later if you want

% % Check if only 'JointMasks' were produced, then assume MaskLab = 'JointMasks'
% Dmsk = dir([acqpath 'Masks*.mat']); Djnt = dir([acqpath 'JointMasks_*.mat']);
% if isempty(Dmsk)&~isempty(Djnt) MaskLab = 'JointMasks'; end

slashes = strfind(acqpath,'\');
Run = acqpath(slashes(end-1)+1:end-1);

% Determine fix parameters, load or get values
if ~isempty(find(fxBools==1))
    if ~exist([UpOneDir(acqpath) 'FixedParsVals.mat']) FLIMGetFixedPars(acqpath); end
    load([UpOneDir(acqpath) 'FixedParsVals.mat']);
end

fitting_method = 1;

try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end

%% Load decays from decays struct. Based on OpeningFunction
% Get the decays for this batch fitting. If 'MaskLab' specified, get only
% those decays. Typically for doing a single generic mask fit for all
% embryos. Otherwise, look for decays for all embryos, of form
% 'decays_Pos3_mask1.mat'
D = dir([acqpath '*decays*' MaskLab '*.mat']);

% Load initial decay
load([acqpath D(1).name]);
decind = ~cellfun('isempty',decay_struct);decind = find(decind); decind = decind(1);
dectime = decay_struct{decind}.time;
dt = dectime(2)-dectime(1);
% Also get laser rep rate from any of the sdts in the folder
Dsdt = subdir([acqpath '\*.sdt']);
[bla,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup(Dsdt(1).name));

% Find IRF sdt file for the day, convert to 'FLIMirf.mat'
% Also get bounds needed for fitting
% Find cal_files folder
if exist([acqpath 'cal_files'])==7 % if cal_files is in acqpath
    calpath = [acqpath '\cal_files\'];
elseif exist([UpOneDir(acqpath) 'cal_files'])==7 % if cal_files in daypath
    calpath  = [UpOneDir(acqpath) 'cal_files\'];
else
    error('Cannot locate cal_files folder. Place in daypath or acqpath, please');
end
Dirfsdt = dir([calpath '*irf*.sdt']);
% Assume only 1 IRF in cal_files
if length(Dirfsdt)>1 warning('Should have only 1 IRF per acquisition! Wassamaddawichu?!'); end
if ~isempty(Dirfsdt)
    IRFstruct = FLIMLoadIRF([calpath Dirfsdt(1).name]);
    irftime = IRFstruct.time;
    [dt_irf,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup([calpath Dirfsdt(1).name])); % Assume the same even if using different IRFs
    % Calc adcratio
    adcratio = length(irftime)/length(dectime);
else
    error('No IRF detected in ''cal_files'' folder');
end

% Figure out how many segments there are
numsegs = size(decay_struct{decind}.decay,2);
if numsegs==3
    segnames = {'-joint','-mito','-cyto'};
else
    segnames{1} = '';
end

% Determine which channels are present
uchans = unique(nameinds(:,4));
if ~isempty(find(strcmp(uchans,'NADH'))) Chind(1) = 1; end
if ~isempty(find(strcmp(uchans,'FAD'))) Chind(2) = 2; end
if ~isempty(find(strcmp(uchans,'UserChan'))) Chind(3) = 3; end
Chind(Chind==0)=[];
ChLabs = {'NADH','FAD','UserChan'};

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

% Initialize the 'mdpars' array (see explanation below where it's filled in).
PosNames = unique(nameinds(:,3));
Tdim = length(unique([nameinds{:,2}]));
Pdim = length(unique(nameinds(:,3)));
Chdim = length(unique(nameinds(:,4)));
Zdim = length(unique([nameinds{:,5}]));
for decnum = 1:size(D,1)
    % If it says 'Joint' anywhere in the title, assume only 1 'joint' mask
    if ~isempty(strfind(D(decnum).name,'Joint'))
        Mskdim = 1;
        Msknums(decnum) = 1;
    else
        mind = strfind(D(decnum).name,'mask');
        Msknums(decnum) = str2num(D(decnum).name(mind+4:end-4));
    end
end
Mdim = max(Msknums);
% if exist([acqpath 'multiD_pars.mat'])==2
%     load([acqpath 'multiD_pars.mat']);
%     if size(mdpars,1)~=N_param+2 error('mdpars size mismatch. Num expos changed?'); end
% else
mdpars = nan(N_param+2,2,Tdim,Pdim,3,Zdim,Mdim,numsegs);
% Note, 'channels' dim always has 3: 1=NADH, 2=FAD, 3=UserChan
% end

% Loop over decay files present
for decnum = 1:size(D,1)
    pind = strfind(D(decnum).name,'Pos');
    dashes = strfind(D(decnum).name,'_');
    PosNum = D(decnum).name(pind+3:dashes(end)-1);
    mind = strfind(D(decnum).name,'mask');
    if length(dashes)>2
        Mind = str2num(D(decnum).name(mind+4:dashes(3)-1));
    else
        Mind = str2num(D(decnum).name(mind+4:end-4));
    end
    if isempty(Mind) Mind = 1; end % Assume single mask if 'mask' not present
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),PosNum)
            continue;
        end
    end
    
    load([acqpath D(decnum).name]);
    
    %     % Sum all decays to have robuts decay for finding fit range
    %     decsum = zeros(256,1);
    %     nonemp = find(~cellfun('isempty',decay_struct));
    %     for dc = 1:length(nonemp) decsum = decsum + sum(decay_struct{nonemp(dc)}.decay,2); end
    %     plot(decsum)
    
    % Separate into channels
    clear chans
    if isempty(find(~cellfun('isempty',decay_struct)))
        continue;
    end
    for k = 1:size(decay_struct,1)
        if ~isempty(decay_struct{k})
            fname = decay_struct{k}.name;
            dashes = strfind(fname,'_');
            chans{k} = fname(dashes(2)+1:dashes(3)-1);
            chnums(k) = find(strcmp(ChLabs,chans{k}));
        end
    end
    
    uchans = unique(chans(~cellfun('isempty',chans)));
    
    % Initialize decay_struct structureson-empty decay_struct
    % elements
    NumOfDecays = length(decay_struct);
    
    names = cell(NumOfDecays,1);
    
    for i = 1:NumOfDecays
        if isempty(decay_struct{i})
            continue;
        end
        
        % Synthetic IRF?
        if exist('IRFstruct')
            decay_struct{i}.irf = IRFstruct.irf;
            decay_struct{i}.IRFdt_ns = dt_irf;
            decay_struct{i}.laserT_ns = LaserT;
            irf = decay_struct{i}.irf;
            DecMax(i) = -1;
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
            sumdec  = sum(decay_struct{i}.decay,2);
            tmp = find(sumdec==max(sumdec)); % offset empirical
            DecMax(i) = tmp(1);
        end
        
        % Get fit and noise regions, if not already done in 'FLIM_decay_from_probmap.m' (updated 2019-05-22).
        if ~isfield(decay_struct{i},'fit_region')
            try
                decay_struct{i}.fit_region = [ceil(IRFstruct.fst/adcratio),floor(IRFstruct.fend/adcratio)];
                decay_struct{i}.noise_region = [ceil(IRFstruct.nst/adcratio),floor(IRFstruct.nend/adcratio)];
            catch
                error('Repeat FLIM_decay_from_probmap to set fit_region, or put an IRF into cal_files')
            end
        end
        
        decay_struct{i}.decay_handle = -99;
        decay_struct{i}.fit_handle = -99;
        decay_struct{i}.residual_handle = -99;
        decay_struct{i}.residual = zeros(length(dectime),numsegs);
        decay_struct{i}.fit_result = zeros(N_param,3,numsegs);
        %    decay_struct{i}.fit_result(1,1) = 5;
        %    decay_struct{i}.fit_result(1,3) = 1;
        decay_struct{i}.Chi_sq(1:3) = 0;
        decay_struct{i}.fitting_method = 0;
        
        decay_struct{i} = orderfields(decay_struct{i});
        
        names{i} = decay_struct{i}.name;
    end
    
    % Proceed with fitting
    if fitting_method == 1   % Do LS fitting
        for j =  1:Bin:NumOfDecays
            % Sum over Bins. Bin photons in the 1st of the binned frames
            % (call that the binned time point). Fit that frame. Then set
            % the other frames in the decay_struct to empty elements, [],
            % so that so that they will be excluded from further plotting and
            % averaging.
            % Note, this was really intended to bin z-positions together.
            % If we want to bin time points, we'll need to write new code.
            decay = zeros(length(dectime),numsegs);
            clear chs;
            irrs = []; irr_stds = []; irr_stderrs = [];
            MaxIter = -1;
            for k = 1:Bin
                if j+k-1>length(decay_struct)
                    continue;
                end
                if isempty(decay_struct{j+k-1})
                    continue;
                end
                name = decay_struct{j+k-1}.filename;
                dashes = strfind(name,'_');
                chs{k} = name(dashes(2)+1:dashes(3)-1);
                decay = decay + decay_struct{j+k-1}.decay;
                irrs = [irrs; decay_struct{j+k-1}.irr];
                irr_stds = [irr_stds; decay_struct{j+k-1}.irr_std];
                irr_stderrs = [irr_stderrs; decay_struct{j+k-1}.irr_std/sqrt(length(find(decay_struct{j+k-1}.selected_pixel)))];
                
                decay_struct{j} = decay_struct{j+k-1}; % copy 'meta data' to binned point
                if k>1 % remove all but the last point in the bin series
                    decay_struct{j+k-1} = [];
                end
            end
            
            if isempty(decay_struct{j})
                continue;
            end
            
            if ~isempty(strfind(decay_struct{j}.filename,'Ruth'))
                continue;
            end
            
            if length(unique(chs(~cellfun('isempty',chs))))>1
                error('Indexing problem. Code is trying to bin NADH and FAD decays together')
            end
            
            % Fill in effective values for binned point
            decay_struct{j}.irr = sum(irrs./irr_stds.^2,1)./sum(1./irr_stds.^2,1);
            decay_struct{j}.irr_std = sqrt(1./sum(1./irr_stds.^2,1));
            decay_struct{j}.irr_stderr = sqrt(1./sum(1./irr_stderrs.^2,1));
            % Note, we're combining irr measures from different planes
            % here (rather than total photons/total pixels from all planes)
            decay_struct{j}.Bin = Bin;
            
            decay_struct{j}.decay = decay;
            time = decay_struct{j}.time;
            fit_start = decay_struct{j}.fit_region(1);
            fit_end = decay_struct{j}.fit_region(2);
            noise_start = decay_struct{j}.noise_region(1);
            noise_end = decay_struct{j}.noise_region(2);
            counts = sum(decay(fit_start:fit_end,:),1);
            decay_struct{j}.photons = counts;
            
            
            %weight on residual
            %weight = (fit_end-fit_start+1)/sqrt(decay(fit_start:fit_end)'*decay(fit_start:fit_end));
            nonzero_decay = decay;
            nonzero_decay(find(decay==0))=1;
            sigy = sqrt(nonzero_decay);
            weight = 1./sigy(fit_start:fit_end,:);
            
            % rough estimate of noise in noise region
            est_noise = mean(decay(noise_start:noise_end,1));
            
            % Update the noise initial estimate, use inputs above for the
            % rest of the parameters.
            p_init = p_init0; p_init(2) = 1-est_noise/max(decay(:,1));
            p_min = p_min0;
            p_max = p_max0;
            
            fit_result = zeros(N_param,3);
            
            % Parameters to fix (shift is main one). Loop over fix options.
            for fx = 1:length(fxBools)
                if fxBools(fx)
                    fit_result(fx,1) = fxvals(fx,chnums(j));
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
            
            % Proceed with fits, but loop over different segments (mito,
            % ctyo, joint), each with their respective decays
            decay_struct{j}.residual = zeros(size(decay));
            decay_struct{j}.fit = zeros(length(dectime(fit_start:fit_end)),numsegs);
            
            for seg = 1:size(decay,2)
                if ~isempty(find(decay(:,seg))) % skip seg fit for any empty segments
                    [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst, converged] = ...
                        lm(@lm_decay_model,p_init,time,decay(:,seg),weight(:,seg),dp,p_min,p_max,[nexpo,counts(seg),fit_start,fit_end,dt_irf,LaserT,DecMax(j)],fit_start,fit_end,irf,MaxIter);
                    
                    
                    %                     % sorting test. Plot fits
                    %                     p_fit = p_fit([1 2 5 4 3]); p_fit(4)=1-p_fit(4);
                    %                     y_hat = lm_decay_model(time,p_fit,[nexpo,counts(seg),fit_start,fit_end,dt_irf,LaserT],irf);
                    %                     y_hat = y_hat(fit_start:fit_end);
                    %                     close all;
                    %                     semilogy(y_hat,'b'); hold on;
                    
                    % Sort params in order of increasing lifetimes.
                    [p_fit,sigma_p] = SrtFLIMparamsByTaus(p_fit,sigma_p);
                    % Note that 'corr' is not sorted. Could do with ij for loops, but eh.
                    
                    disp([Run '\' D(decnum).name(8:end-4) '-Fr' num2str(nameinds{j,6}) '-' chans{j} segnames{seg} ': ' num2str(size(cvg_hst,1)) ' iterations'])
                    
                    y_hat = lm_decay_model(time,p_fit,[nexpo,counts(seg),fit_start,fit_end,dt_irf,LaserT,DecMax],irf);
                    y_hat = y_hat(fit_start:fit_end);
                    
                    decay_struct{j}.fit(:,seg) = y_hat;
                    
                    y_dat = decay(:,seg);
                    y_dat = y_dat(fit_start:fit_end);
                    %                     close all; semilogy(y_hat); hold on; plot(y_dat,'r') % To check fit
                    
                    weighted_residual = weight(:,seg).*(y_dat(:)-y_hat(:));
                    decay_struct{j}.residual(fit_start:fit_end,seg) = weighted_residual;
                    
                    fit_result(1:N_param,1) = real(p_fit);
                    fit_result(1:N_param,2) = real(sigma_p);
                    
                    decay_struct{j}.fit_result(:,:,seg) = fit_result;
                    decay_struct{j}.Chi_sq(seg) = Chi_sq;
                    %                 decay_struct{j}.cvg_hst(seg) = cvg_hst;
                    
                    decay_struct{j}.fit_region = [fit_start,fit_end];
                    decay_struct{j}.noise_region = [noise_start,noise_end];
                    decay_struct{j}.fitting_method = fitting_method;
                    decay_struct{j}.nexpo = nexpo;
                    decay_struct{j}.converged = converged;
                    
                    %                     % Test plots
                    %                     close all;figure('position',[100 100 450 450])
                    %                     sp1=subplot(2,1,1);
                    %                     semilogy(time(fit_start:fit_end),decay(fit_start:fit_end,seg)); hold on;
                    %                     plot(time(fit_start:fit_end),y_hat)
                    %                     set(gca,'xticklabels',[]); axis tight;
                    %                     sp2=subplot(2,1,2);
                    %                     plot(time(fit_start:fit_end),weighted_residual); hold on;
                    %                     plot([time(fit_start) time(fit_end)],[0 0],'k-.')
                    %                     sp1.Position = [0.13 0.32 0.775 0.605];
                    %                     sp2.Position = [0.13 0.14 0.775 0.167];
                    %                     xlabel('time (ns)'); axis tight; ylim([-5 5])
                    
                    % Build 'mdpars', a multi-dimensional (md) array that
                    % contains all the parameters in one place. Dimensions
                    % start with parameter indices, then dims 3:end are in
                    % the same order as nameinds:
                    % [Param#, mean/stderr(1,2), time point, position, channel, z-position, embryo number, segment]
                    % So 8 in total, can add more later if we need to.
                    fname = decay_struct{j}.name;
                    [Tind,Chind,Zind] = MultiDindsFromSdtName(fname);
                    if Zdim==Bin Zind = 1; end
                    Pind = find(strcmp(PosNames,PosNum));
                    irr = decay_struct{j}.irr;
                    irr_stderr = decay_struct{j}.irr_stderr;
                    if j==20
                        1;
                    end
                    
                    mdpars(:,:,Tind,Pind,Chind,Zind,Mind,seg)=[fit_result(:,1:2);[irr(seg) irr_stderr(seg)];[decay_struct{j}.timestp decay_struct{j}.timestp]];
                    % Bad frames remain filled in with NaNs
                end
            end
            decay_struct{j} = orderfields(decay_struct{j});
        end
    elseif fitting_method == 3
        %         % Should work, but will have to play with fitting bounds or else it
        %         % takes forever. Ask Tae, and stick with LSq for now
        %         %     matlabpool;    %for parallel computing
        %         prior = 1;
        %         for j = 1:NumOfDecays
        %
        %             time = decay_struct{j}.time;
        %
        %             %parameters needed for fitting
        %             dt = time(2)-time(1);
        %             fit_start = round(fit_start_x/dt);
        %             fit_end = round(fit_end_x/dt);
        %
        %
        %             % data to fit
        %             decay = decay_struct{j}.decay;
        %             counts = sum(decay(fit_start:fit_end));
        %
        %
        %             %weight on residual
        %             %weight = (fit_end-fit_start+1)/sqrt(decay(fit_start:fit_end)'*decay(fit_start:fit_end));
        %             nonzero_decay = decay;
        %             nonzero_decay(decay==0)=1;
        %             sigy = sqrt(nonzero_decay);
        %             weight = 1./sigy(fit_start:fit_end);
        %
        %             % initial guess, lower and upper bounds of parameters
        %             % P(t|param) (not normalized)
        %             % if one expo:
        %             % P = param(2)*exp(-t/param(3))+(1-param(2))
        %             % param(2): fractional amp of expo decay
        %             % param(3): lifetime of the first decay
        %             %
        %             % if two expo:
        %             % P = param(2)*param(4)*exp(-t/param(3))+param(2)*(1-param(4))*exp(-t/(param(5)*param(3)))+(1-param(2));
        %             % param(4) : fraction of the first decay
        %             % param(5) : 1-FRET efficiency (short/long lifetime)
        %             %
        %             % param(1) : shift of decay model from IRF (usually ranges from -10 to 10)
        %             if nexpo == 1
        %                 p_min = [-15,0.96,2]';
        %                 p_max = [-5,1,3]';
        %                 dp = [1,0.00025,0.01]';
        %             elseif nexpo == 2
        %                 p_min = [-15,0.96,2.1,0.3,0.01]';
        %                 p_max = [-5,1,3,1,0.3]';
        %                 dp = [1,0.0025,0.01,0.01,0.01]';
        %             else
        %                 errordlg('not available yet')
        %                 return
        %             end
        %
        %             % Parameters to fix
        %             fit_result = decay_struct{j}.fit_result;
        %             fixed = fit_result(1:(2*nexpo+1),3);
        %             p_min(fixed==1) = fit_result(fixed==1,1);
        %             p_max(fixed==1) = fit_result(fixed==1,1);
        %
        %
        % %             tic;
        %             [p_fit,sigma_p,p_vec,post,marg_post] = bayes_fit(time,decay,dp,p_min,p_max,nexpo,prior,fit_start,fit_end,0,1);
        % %             toc;
        %
        %             y_hat = lm_decay_model(time,p_fit,[nexpo,counts,fit_start,fit_end]);
        %             y_hat = y_hat(fit_start:fit_end);
        %
        %             decay_struct{j}.fit = y_hat;
        %
        %             y_dat = decay;
        %             y_dat = y_dat(fit_start:fit_end);
        %
        %             weighted_residual = weight.*(y_dat-y_hat);
        %
        %             decay_struct{j}.residual = zeros(size(decay));
        %             decay_struct{j}.residual(fit_start:fit_end) = weighted_residual;
        %             decay_struct{j}.residual_handle = plot(time(fit_start:fit_end),weighted_residual);
        %
        %             Chi_sq = sum(weighted_residual.^2)/(fit_end-fit_start-2*nexpo-1+sum(fixed));
        %
        %             fit_result(1:5,1:2) = zeros(5,2);
        %             fit_result(1:(nexpo*2+1),1) = real(p_fit);
        %             fit_result(1:(nexpo*2+1),2) = real(sigma_p);
        %
        %             decay_struct{j}.fit_result = fit_result;
        %             decay_struct{j}.Chi_sq = Chi_sq;
        %             decay_struct{j}.marg_post = marg_post;
        %             decay_struct{j}.p_vec = p_vec;
        %             decay_struct{j}.post = post;
        %
        %             decay_struct{j}.fit_region = [fit_start,fit_end];
        %             decay_struct{j}.fitting_method = fitting_method;
        %             decay_struct{j}.nexpo = nexpo;
        %             decay_struct{j}.prior = prior;
        %
        %             decay_struct{j} = orderfields(decay_struct{j});
        %         end
        %     matlabpool close
    end
    
    
    % Save fits as a session that can be loaded into the GUI later for viewing
    % Based on 'SaveSession_pushbutton_Callback'
    %     if fxBools(1)
    %         fitname = ['fits' D(decnum).name(7:end-4) '_fxshft.mat'];
    %     else
    fitname = ['fits' D(decnum).name(7:end-4) '.mat'];
    %     end
    file = [acqpath fitname];
    decays_fits_struct = decay_struct;
    save(file,'decays_fits_struct')
    
    
    % Get photon counts to put in spreadsheet for quick check.
    tmp = decays_fits_struct(~cellfun('isempty',decays_fits_struct));
    dashes = fitname;
    MasLab = fitname(dashes(2)+1:strfind(fitname,'.mat')-1);
    spr{1,1} = 't-point';
    
    spr{1,2*(decnum-1)+2} = ['N_P' num2str(PosNum) 'm' MasLab];
    spr{1,2*(decnum-1)+3} = ['F_P' num2str(PosNum) 'm' MasLab];
    spr{1,2*(decnum-1)+4} = ['U_P' num2str(PosNum) 'm' MasLab];
    for i = 1:length(tmp)
        dashes = strfind(tmp{i}.filename,'_');
        t = str2num(tmp{i}.filename(dashes(1)+1:dashes(2)-1))+1;
        spr{t+1,1}=t;
        spr{t+1,chnums(i)+1} = sum(tmp{i}.decay(:,1));
    end
    
    %     if nexpo==2 FLIMBatchParamPlot(file); end % Replaced with 'FLIMAcqParamPlot' below
    clear decay_struct decays_fits_struct
end

xlswrite([acqpath 'PhotCounts.xls'],spr);

if strcmp(MaskLab,'mask')
    save([acqpath 'multiD_pars.mat'],'mdpars');
else
    save([acqpath 'multiD_pars_' MaskLab '.mat'],'mdpars');
end

% MEAN
% ts = unique(nanmean(nanmean(nanmean(mdpars(7,1,:,:,:,:,1,1),4),5),6)); ts(isnan(ts))=[];
% if TplotBool&length(ts)>1 FLIMAcqParamPlot(acqpath); end
FLIMAcqParamPlot(acqpath);
