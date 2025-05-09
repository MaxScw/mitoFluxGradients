function LaserOPplot(acqs,savedest)
% Read in AcqSdtCorr.txt files, extract logged laser output power values,
% and plot them. Also save 'laserOPs.mat' if you want to later scale
% intensities by power^2 (for 2-photon).
% INPUTS:
% -acqs: a cell of acquisition paths.
% -savedest: path and filename for output

% clear all; close all;
% daypath = 'Z:\Lab\Marta\Flow Chamber Tests\2019-05-07 Rotenone w oxamate pipette\';
% acqs{1} = [daypath 's1_a1\'];
% acqs{2} = [daypath 's1_a1b\'];
% acqs{3} = [daypath 's1_a1c\'];
% savedest = [daypath 'LaserOPs'];

AcqSdtCorr = [];
for i = 1:length(acqs)
    AcqSdtCorr = [AcqSdtCorr; ReadTxtSprd2Cell([acqs{i} 'AcqSdtCorr.txt'])];
end

% Find channels
chans = AcqSdtCorr(:,4);
uchans = unique(chans); 
uchans = flip(uchans); % So NADH typically comes first
minT = min(datenum(AcqSdtCorr(:,6)));
cols = [[0 0 1];[0 1 0];[1 0 0]];
% Plot pow for each channel
figure('position',[200 200 450*length(uchans) 400])
for ch = 1:length(uchans)
    clear t pow
    powcl = AcqSdtCorr(strcmp(AcqSdtCorr(:,4),uchans(ch)),end);
    t = datenum(AcqSdtCorr(strcmp(AcqSdtCorr(:,4),uchans(ch)),6));
    t = (t-minT)*86400/60; % min
    for i = 1:length(powcl)
        pow(i) = str2num(powcl{i}(end-7:end-1));
    end
    [a,b] = sort(t); t=t(b); pow=pow(b);
    subplot(1,length(uchans),ch)
    plot(t,pow,'color',cols(:,ch));
    xlabel('time (min)'); ylabel([uchans(ch) 'Laser output pow (W)']); %legend(uchans(ch),'location','best')
    Ts{ch} = t; Pows{ch} = pow;
end
if exist('savedest')&savedest~=-1
    saveas(gcf,[savedest '.fig']);
    save([savedest '.mat']);
end