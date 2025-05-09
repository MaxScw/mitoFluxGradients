function FLIM_batch_fitting_weka(path,poss,fxshftBool,fxBgBool,IRFname,Bin,MaskLab,nexpo)
% Simple. Search path for .mat files with 'decays' in them, perform fits
% for all contained decays and save as 'fits_...'
% Take guts from 'FittingGUI_ver3_4'
% INPUTS
% -fxshftBool: fix shifts option, using 'shifts.mat' in root directory,
%  obtained from 'FLIMGetFixedShifts.m'.
% -IRFname: file name of IRF sdt file stored in the root directory. If
%  blank, function just looks for an sdt with 'irf' in it
% -poss: optional, fit for only one position.
% -Bin: Bins together every 'Bin' frames into  a single decay and then
%   fits. Useful if photon count is low.


% Versions:
% 2018-08-11: Major rewrite - implement 3-exp fitting option
% 2018-08-03: Not using the Weka segmentation anymore, but I had to update
% this version to correct the TAC problems in lm_decay_model.m Used to
% reanayze CLPP data for Emre collaboration, after which we switched to
% Ilastik probability maps approach.
% 2018-01-18: incorporated new 'cal_files' organization
% 2015-05-29: Add binning capability.
%  NOTE: ONLY USEFUL FOR BINNING Z RIGHT NOW. IF Bin>#Z, IT WILL AVERAGE
%   CONSECUTIVE NADH AND FAD FRAMES.
%   If I ever want to come back to this, start with this code:
%   zbin = 3-1; tbin = 2;
%   uManPosNum = 1;
%   Ninds = find(([nameinds{:,2}]==uManPosNum)&strcmp(nameinds(:,4),'NADH')'&([nameinds{:,7}]>-1));
%   Ninds_zbin = find(([nameinds{:,2}]==uManPosNum)&([nameinds{:,5}]==zbin)&strcmp(nameinds(:,4),'NADH')'&([nameinds{:,7}]>-1));
%   length(find(Ninds))
%   length(find(Ninds_zbin))
% 2015-04-17: Adapt to include two separate IRFs for the NADH and FAD fits,
%   respectively. Find IRFs with appropriate names.
% 2015-02-17: Set the fit regions to be automatically detected. Borrow from
%   anisotropy code - "SyncPolChans4GUI"
% 2014-11-05: Include additional fixed-shift option. Calculates shift from
%   first 5 frames binned for each channel, then uses that fixed shift for
%   all frames. Also remove Bayes section for now. It's just confusing. Add
%   back in later if I ever need to. fxshftBool is 1 if you want to fix the
%   shifts. Must previously run 'FLIMGetFixedShifts', which calculates the
%   shifts and saves in root directory.

% clear all;
% path = 'C:\Users\Tim\Documents\Academic - Research\Data\REWRITE2EXP\s1_a1_Day1_test';
% Bin = 3;
% fxshftBool=1;
% nexpo =2;
% fxBgBool = 1;

if path(end)~='\' path = [path '\']; end;
if ~exist('fxshftBool')|fxshftBool==-1 fxshftBool=0; end
if ~exist('fxBgBool')|fxBgBool==-1 fxBgBool=0; end
if ~exist('Bin')|Bin==-1 Bin=1; end
if ~exist('MaskLab')|MaskLab==-1 MaskLab='mask'; end
if ~exist('nexpo')|nexpo==-1 nexpo=2; end

% nexpo = 2; % Change later if you want

slashes = strfind(path,'\');
Run = path(slashes(end-1)+1:end-1);

if fxshftBool|fxBgBool load([UpOneDir(path) 'ShiftsBG.mat']); end
fitting_method = 1;

try     load([path 'multiD_indices.mat']); catch     load([path 'name_indexes.mat']); end

%% Load decays from decays struct. Based on OpeningFunction
% Get the decays for this batch fitting. If 'MaskLab' specified, get only
% those decays. Typically for doing a single generic mask fit for all
% embryos. Otherwise, look for decays for all embryos, of form
% 'decays_Pos3_mask1.mat'
D = dir([path '*decays*' MaskLab '*.mat']);

% Find IRF sdt file for the day, convert to 'FLIMirf.mat'
% Also get bounds needed for fitting
% Check for IllProfCal.mat file
if exist([path 'cal_files'])==7 % if new cal_files folder exists
    CalPath = [path 'cal_files\'];
else
    CalPath = [UpOneDir(path) 'DailyFiles\'];
end
Dirfsdt = dir([CalPath '*irf_*.sdt']);
if exist('IRFname')&IRFname~=-1
    IRFstruct = FLIMLoadIRF([CalPath IRFname]);
else
    for i = 1:length(Dirfsdt)
        % NADH IRF taken at 750
        if ~isempty(strfind(Dirfsdt(i).name,'750nm'))
            NIRFstruct = FLIMLoadIRF([CalPath Dirfsdt(i).name]);
            irftime = NIRFstruct.time;
            % Calculate time bin size of the irf, and get laser rep rate (period, actually)
        end
        % FAD IRF taken at 845
        if ~isempty(strfind(Dirfsdt(i).name,'890nm'))|~isempty(strfind(Dirfsdt(i).name,'845nm'))
            FIRFstruct = FLIMLoadIRF([CalPath Dirfsdt(i).name]);
            irftime = FIRFstruct.time;
        end
    end
end


% For new setup, the 750nm IRF is messy and not that offset from the 845nm
% IRF. Therefore, using the 845nm IRF for all fits gives more stable FLIM
% results
NIRFstruct = FIRFstruct;
[dt_irf,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup([CalPath Dirfsdt(i).name])); % Assume the same even if using different IRFs

% Calc adcratio
load([path D(1).name]);
decind = ~cellfun('isempty',decay_struct);decind = find(decind); decind = decind(1);
dectime = decay_struct{decind}.time;
dt = dectime(2)-dectime(1);
adcratio = length(irftime)/length(dectime);

% Figure out how many segments there are
numsegs = size(decay_struct{decind}.decay,2);
if numsegs==3
    segnames = {'-joint','-mito','-cyto'};
else
    segnames{1} = '';
end

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

N_param = 2*nexpo+1;
if nexpo == 1
    %         [shf,     A,    tau1]
    p_init0 = [0,       1,       1]';
    p_min0 =  [-10,   0.7,    0.05]';
    p_max0 =  [10,      1,       5]';
elseif nexpo == 2
    %         [shf,     A,    tau1,      f2,    tau2]
    p_init0 = [  0,     1,     0.2,     0.3,       3]';
    p_min0 =  [-200,  0.7,    0.01,    0.01,    0.05]';
    p_max0 =  [200,     1,     0.9,       1,       5]';
    
%     % Old stuff.
%     p_init0 = [0,1,3,0.3,0.2]';
%     p_min0 = [-200,0.7,0.05,0.01,0.01]';
%     p_max0 = [200,1,5,1,0.9]';
elseif nexpo == 3
    %         [ shf,    A,    tau1,      f2,    tau2,     f3,    tau3]
    p_init0 = [   0,    1,      .1,     0.3,       1,     .3,       3]';
    p_min0 =  [-200,  0.7,    0.01,    0.01,    0.1,   0.01,       1]';
    p_max0 =  [ 200,    1,      2,       1,       4,      1,       6]';
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
%     if ~isempty(strfind(D(decnum).name,'SingleMasks'))
%         Mskdim = 1;
    if ~isempty(strfind(D(decnum).name,'mask'))
        mind = strfind(D(decnum).name,'mask');
        Msknums(decnum) = str2num(D(decnum).name(mind+4:end-4));
    else
        Mskdim = 1; % Assume custom label and single masks.
        Msknums(decnum) = 1;
    end
end
Mdim = max(Msknums);
% if exist([path 'multiD_pars.mat'])==2
%     load([path 'multiD_pars.mat']);
%     if size(mdpars,1)~=N_param+2 error('mdpars size mismatch. Num expos changed?'); end
% else
    mdpars = nan(N_param+2,2,Tdim,Pdim,Chdim,Zdim,Mdim,numsegs);
% end

% Loop over decay files present
for decnum = 1:size(D,1)
    pind = strfind(D(decnum).name,'Pos');
    dashes = strfind(D(decnum).name,'_');
    uManPos = D(decnum).name(pind+3:dashes(end)-1);
    if ~isempty(strfind(D(decnum).name,'mask'))
        mind = strfind(D(decnum).name,'mask');
        if length(dashes)>2
            Mind = str2num(D(decnum).name(mind+4:dashes(3)-1));
        else
            Mind = str2num(D(decnum).name(mind+4:end-4));
        end
    else
        Mind = 1;
    end
    %     strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
    %     uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),uManPos)
            continue;
        end
    end
    
    load([path D(decnum).name]);
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
        end
    end
    
    uchans = unique(chans(~cellfun('isempty',chans)));
    
    % Proceed with fitting and saving for non-empty decay_struct
    % elements
    NumOfDecays = length(decay_struct);
    
    names = cell(NumOfDecays,1);
    
    for i = 1:NumOfDecays
        if isempty(decay_struct{i})
            continue;
        end
        
        % Figure out which irf to use
        if strcmp(chans{i},'NADH')
            decay_struct{i}.irf = NIRFstruct.irf;
            decay_struct{i}.fit_region = [ceil(NIRFstruct.fst/adcratio),floor(NIRFstruct.fend/adcratio)];
            decay_struct{i}.noise_region = [ceil(NIRFstruct.nst/adcratio),floor(NIRFstruct.nend/adcratio)];
        elseif strcmp(chans{i},'FAD')
            decay_struct{i}.irf = FIRFstruct.irf;
            decay_struct{i}.fit_region = [ceil(FIRFstruct.fst/adcratio),floor(FIRFstruct.fend/adcratio)];
            decay_struct{i}.noise_region = [ceil(FIRFstruct.nst/adcratio),floor(FIRFstruct.nend/adcratio)];
        elseif strcmp(chans{i},'Ruth')
            continue;
        elseif strcmp(chans{i},'UserChan') % if custom user channel, e.g. 900nm illumination, intensity is typically only param of interest.
            decay_struct{i}.irf = FIRFstruct.irf;
            decay_struct{i}.fit_region = [ceil(FIRFstruct.fst/adcratio),floor(FIRFstruct.fend/adcratio)];
            decay_struct{i}.noise_region = [ceil(FIRFstruct.nst/adcratio),floor(FIRFstruct.nend/adcratio)];
        else
            error('Channel problem. Not FAD or NADH')
        end
        decay_struct{i}.decay_handle = -99;
        decay_struct{i}.IRFdt_ns = dt_irf;
        decay_struct{i}.laserT_ns = LaserT;
        decay_struct{i}.fit_handle = -99;
        decay_struct{i}.residual_handle = -99;
        decay_struct{i}.residual = [];
        decay_struct{i}.fit_result = zeros(N_param,3);
        %    decay_struct{i}.fit_result(1,1) = 5;
        %    decay_struct{i}.fit_result(1,3) = 1;
        decay_struct{i}.Chi_sq = 0;
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
            
            %             if ~isempty(strfind(decay_struct{j}.filename,'UserChan'))
            %                 continue;
            %             end
            
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
            irf = decay_struct{j}.irf;
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

            % Parameters to fix (shift is main one)
            if fxshftBool
                if strfind(decay_struct{j}.name,'NADH')
                    fit_result(1,1) = shifts(1);
                elseif strfind(decay_struct{j}.name,'FAD')
                    fit_result(1,1) = shifts(2);
                elseif strfind(decay_struct{j}.name,'UserChan')
                    fit_result(1,1) = shifts(2);
                else
                    error('Channel error')
                end
                fit_result(1,3) = 1;
            end
            if fxBgBool
                if strfind(decay_struct{j}.name,'NADH')
                    fit_result(2,1) = As(1);
                elseif strfind(decay_struct{j}.name,'FAD')
                    fit_result(2,1) = As(2);
                elseif strfind(decay_struct{j}.name,'UserChan')
                    fit_result(2,1) = As(2);
                else
                    error('Channel error')
                end
                fit_result(2,3) = 1; 
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
            
            % Proceed with fits, but loop over different segments (mito,
            % ctyo, joint), each with their respective decays
            decay_struct{j}.residual = zeros(size(decay));
            decay_struct{j}.fit = zeros(length(dectime(fit_start:fit_end)),numsegs);
            for seg = 1:size(decay,2)
                if ~isempty(find(decay(:,seg))) % skip seg fit for any empty segments
                    [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst, converged] = ...
                        lm(@lm_decay_model,p_init,time,decay(:,seg),weight(:,seg),dp,p_min,p_max,[nexpo,counts(seg),fit_start,fit_end,dt_irf,LaserT],fit_start,fit_end,irf,MaxIter);
                    
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
                    
                    y_hat = lm_decay_model(time,p_fit,[nexpo,counts(seg),fit_start,fit_end,dt_irf,LaserT],irf);
                    y_hat = y_hat(fit_start:fit_end);
                    %                     semilogy(y_hat,'r.')
                    
                    decay_struct{j}.fit(:,seg) = y_hat;
                    
                    y_dat = decay(:,seg);
                    y_dat = y_dat(fit_start:fit_end);
                    
                    weighted_residual = weight(:,seg).*(y_dat(:)-y_hat(:));
                    decay_struct{j}.residual(fit_start:fit_end,seg) = weighted_residual;
                    
                    fit_result(1:(nexpo*2+1),1) = real(p_fit);
                    fit_result(1:(nexpo*2+1),2) = real(sigma_p);
                    
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
                    % the same order as nameinds::
                    % [Param#, mean/std(1,2), time point, position, channel, z-position, embryo number, segment]
                    %  Params: First elements are the output of the fit.
                    %          Then concat irr and the timestamp onto the
                    %          end (std err of irr and timmestamp again
                    %          onto the end of the sigma array).
                    %          e.g. 2-exp: [shft,A,tau1(short),frac,tau2,irr,timestamp]
                    % So 8 in total, can add more later if we need to.
                    fname = decay_struct{j}.name;
                    [Tind,Chind,Zind] = MultiDindsFromSdtName(fname);
                    if Zdim==Bin Zind = 1; end
                    Pind = find(strcmp(PosNames,uManPos));
                    irr = decay_struct{j}.irr;
                    irr_stderr = decay_struct{j}.irr_stderr;
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
%             %             tic;
%             [p_fit,sigma_p,p_vec,post,marg_post] = bayes_fit(time,decay,dp,p_min,p_max,nexpo,prior,fit_start,fit_end,0,1);
%             %             toc;
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
%     if fxshftBool
%         fitname = ['fits' D(decnum).name(7:end-4) '_fxshft.mat'];
%     else
        fitname = ['fits' D(decnum).name(7:end-4) '.mat'];
%     end
    file = [path fitname];
    decays_fits_struct = decay_struct;
    save(file,'decays_fits_struct')
    
    
    % Get photon counts to put in spreadsheet for quick check.
    tmp = decays_fits_struct(~cellfun('isempty',decays_fits_struct));
    dashes = fitname;
    MasLab = fitname(dashes(2)+1:strfind(fitname,'.mat')-1);
    spr{1,1} = 't-point';
    spr{1,2*(decnum-1)+2} = ['N_P' num2str(uManPos) 'm' MasLab];
    spr{1,2*(decnum-1)+3} = ['F_P' num2str(uManPos) 'm' MasLab];
    for i = 1:length(tmp)
        dashes = strfind(tmp{i}.filename,'_');
        t = str2num(tmp{i}.filename(dashes(1)+1:dashes(2)-1))+1;
        spr{t+1,1}=t;
        if ~isempty(strfind(tmp{i}.filename,'NADH'))
            spr{t+1,2*(decnum-1)+2} = sum(tmp{i}.decay);
        elseif ~isempty(strfind(tmp{i}.filename,'FAD'))
            spr{t+1,2*(decnum-1)+3} = sum(tmp{i}.decay);
        else
            error('Not NADH or FAD!')
        end
    end
    
    if nexpo==2 FLIMBatchParamPlot(file); end
    clear decay_struct decays_fits_struct
    
end
xlswrite([path 'PhotCounts.xls'],spr);
save([path 'multiD_pars.mat'],'mdpars');

 