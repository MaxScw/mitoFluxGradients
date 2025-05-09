function FLIMGetFixedPars(acqpath,nexpo) %,fit_start,fit_end,noise_region_from,noise_region_to,p_init,p_min,p_max,irf);
% % Expanded from 'FLIMGetFixedShifts.m' to also calculate average A=(1-BG)
% parameter.
% Finds the average shifts for a given data set. Intended to do on one
% sample, then plug in those NADH and FAD shifts into FLIM_batch_fitting
% for all the samples taken on that day. 
% acqpath - points to a subfolder with a representative sample

% % TEST in script mode:
% clear all;
% acqpath = 'C:\Dropbox\data\s1_a1';

if acqpath(end)~='\' acqpath = [acqpath '\']; end
if ~exist('nexpo')|nexpo==-1 nexpo=2; end
D = dir([acqpath '*decays*.mat']);
% load([acqpath D(1).name]);
% decind = ~cellfun('isempty',decay_struct); decind = find(decind); decind = decind(1);
% timebins = length(decay_struct{decind}.decay);
% time = decay_struct{decind}.time;

% Load initial decay
load([acqpath D(1).name]);
decind = ~cellfun('isempty',decay_struct);decind = find(decind); decind = decind(1);
dectime = decay_struct{decind}.time;
timebins = length(decay_struct{decind}.decay);
dt = dectime(2)-dectime(1);
time = decay_struct{decind}.time;

% Also get laser rep rate from any of the sdts in the folder
Dsdt = subdir([acqpath '\*.sdt']);
[bla,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup(Dsdt(1).name));


%%%%%% Fitting stuff
% Load IRF in root acqpath
% Find cal_files folder
if exist([acqpath 'cal_files'])==7 % if cal_files is in acqpath
    calpath = [acqpath 'cal_files\'];
elseif exist([UpOneDir(acqpath) 'cal_files'])==7 % if cal_files in daypath
    calpath  = [UpOneDir(acqpath) 'cal_files\'];
else
    error('Cannot locate cal_files folder. Place in daypath or acqpath, please');
end
Dirfsdt = dir([calpath '*irf*.sdt']);
% Assume only 1 IRF in cal_files
if length(Dirfsdt)>1 error('Should have only 1 IRF per acquisition! Wassamaddawichu?!'); end
if ~isempty(Dirfsdt)
    IRFstruct = FLIMLoadIRF([calpath Dirfsdt(1).name]);
    irftime = IRFstruct.time;
    [dt_irf,LaserT] = GetLaserT_dt_FromSdt(bh_readsetup([calpath Dirfsdt(1).name])); % Assume the same even if using different IRFs
    % Calc adcratio
    adcratio = length(irftime)/length(dectime);
    DecMax = -1;
else
    % If no IRF present, lm_decay creates a synthetic gaussian IRF,
    % but it helps to do an initial centering by finding the knee
    sumdec  = sum(decay_struct{decind}.decay,2);
    tmp = find(sumdec==max(sumdec)); % offset empirical
    DecMax = tmp(1);
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
    p_min0 =  [-200,  0.7,    0.01,    0.01,    0.05]';
    p_max0 =  [200,     1,       1,       1,      10]';
elseif nexpo == 3
    %         [ shf,    A,    tau1,      f2,    tau2,     f3,    tau3]
    p_init0 = [   0,    1,      .2,     0.3,     0.5,     .3,       3]';
    p_min0 =  [-200,  0.7,    0.05,    0.01,     1,   0.01,       1]';
    p_max0 =  [ 200,    1,      .1,       1,       5,      1,      10]';
else
    errordlg('not available yet')
end

% Loop over sdt's. For each decay_struct file, bin all the elements
% together, to get one big decay curve. Each channel is binned separately,
% though.
for dec = 1:size(D,1)
    load([acqpath D(dec).name]);
    % One 'decays array' with three channels present in 3 columns: 1=NADH, 2=FAD, 3=User
    decays = zeros(timebins,3);
    ind = find(~cellfun('isempty',decay_struct))';
    fit_start = decay_struct{ind(1)}.fit_region(1);
    fit_end = decay_struct{ind(1)}.fit_region(2);
    noise_start = decay_struct{ind(1)}.noise_region(1);
    noise_end = decay_struct{ind(1)}.noise_region(2);
    
    % Synthetic IRF?
    if exist('IRFstruct')
        irf = IRFstruct.irf;
    else
        irf = -1; dt_irf = -1;
        % For synthetic IRF, fill in the last fit parameter as gaussian width
        % Initial value for IRF width measured empirically from a
        % typical IRF. ~90 bins on a 4096 bin array
        p_init0(N_param) = 100;
        p_min0(N_param) = 1;
        p_max0(N_param) = 300;
    end

    % Loop over all decay elements and bin photons into one decay for each
    % channel
    for j = ind
        if ~isempty(strfind(decay_struct{j}.name,'NADH'))
            decays(:,1) = decays(:,1) + decay_struct{j}.decay(:,1); 
            % only use 'joint' segment (mito + cyto together)
        elseif ~isempty(strfind(decay_struct{j}.name,'FAD'))
            decays(:,2) = decays(:,2) + decay_struct{j}.decay(:,1);
        elseif ~isempty(strfind(decay_struct{j}.name,'UserChan'))
            decays(:,3) = decays(:,3) + decay_struct{j}.decay(:,1);
        end
    end
    
    % Loop over channels and fit binned decays to get pars
    for ch = 1:3
        if ~isempty(find(decays(:,ch))) 
            
            nonzero_decay = decays(:,ch);
            nonzero_decay(decays(:,ch)==0)=1;
            sigy = sqrt(nonzero_decay);
            weight = 1./sigy(fit_start:fit_end);
            est_noise = mean(decays(noise_start:noise_end,1));
            p_init0(2) = 1-est_noise/max(decays(:,ch));
            dp = ones(N_param,1).*0.001; dp(1) = 1;
            if ~exist('IRFstruct') dp(end) = 1; end
            [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst, converged,y_hat] = ...
                lm(@lm_decay_model,p_init0,time,decays(:,ch),weight,dp,p_min0,p_max0,[nexpo,sum(decays(fit_start:fit_end,ch)),fit_start,fit_end,dt_irf,LaserT,DecMax],fit_start,fit_end,irf);
%             shifts(dec,ch) = p_fit(1);
%             As(dec,ch) = p_fit(2);
            fxvals(:,dec,ch) = p_fit;
%             if exist('IRFstruct')
%                 SynthIrfW(dec,ch) = -1;
%             else
%                 SynthIrfW(dec,ch) = p_fit(end);
%             end
            
%             % Test fit:
%             figure;
%             y_hat = lm_decay_model(time,p_fit,[nexpo,sum(decays(fit_start:fit_end,ch)),fit_start,fit_end,dt_irf,LaserT,DecMax],irf); y_hat = y_hat(fit_start:fit_end);
%             y_dat = decays(:,ch); y_dat = y_dat(fit_start:fit_end);
%             semilogy(time(fit_start:fit_end),y_dat,'b',time(fit_start:fit_end),y_hat,'k');
            1;
        end
    end
end

% Now we have parameter arrays with one value per decay file (corresponding
% to one mask, all frames). Take the averages for fix values.
% shifts = mean(shifts,1);
% As = mean(As,1);
% SynthIrfW = mean(SynthIrfW,1);
fxvals = squeeze(mean(fxvals,2));
save([UpOneDir(acqpath) 'FixedParsVals.mat'],'fxvals')
