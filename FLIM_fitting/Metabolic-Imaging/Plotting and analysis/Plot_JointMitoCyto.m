% Script to take a fits files and make comparative plots of all 3 segments,
% 'mito, 'cyto', and 'joint'. Plots are saved to the same location as the fits 

clear all
% USER INPUT
file = 'C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\2016-11-10 1-cell\s1_a1_O2drop\fits_Pos0_SingleMasks.mat';


datpath = UpOneDir(file);
load(file);


% Optional range of sdt #'s. Comment out if unwanted.
% rnge = 9:50;

% Titles and file names
titles = {'Intensities','Fractions Bound','Short Lifetimes','Long Lifetimes'};
ylabs = {'Int','fbound','\tau_1 (ns)','\tau_2 (ns)'};
savelabs = {'Ints','bounds','tau1s','tau2s'};
suffix = '';

L = length(decays_fits_struct);
% Get NADH and FAD indices
Nind = zeros(L,1); Find = zeros(L,1);
for i = 1:L
    if ~isempty(strfind(decays_fits_struct{i}.name,'NADH')) Nind(i)=1; end
    if ~isempty(strfind(decays_fits_struct{i}.name,'FAD')) Find(i)=1; end
end

% Build the parameters matrix
for i = 1:length(decays_fits_struct) % loop over cell elements (std's)
    if ~isempty(decays_fits_struct{i})
        restab = decays_fits_struct{i}.fit_result;
        irrs = decays_fits_struct{i}.irr';
        fbounds = squeeze(restab(4,1,:));
        tau1s = squeeze(restab(3,1,:)).*squeeze(restab(5,1,:));
        tau2s = squeeze(restab(3,1,:));
        
        % params has all the data. Indices: (param,seg,t)
        params(:,:,i) = [irrs fbounds tau1s tau2s]';
        scs(i) = decays_fits_struct{i}.numscans;
        ts(i) =  decays_fits_struct{i}.timestp;
    end
end

% Get time array from abs time stamps
tsnon0 = ts(ts>0);
ts0 = (ts-min(tsnon0))*86400/60; % minutes

%%
close all
cols = {'r','g','b'};
chans = {'NADH','FAD'};
for ch = 1:2
    
    for p = 1:length(titles)
        h(p) = figure;
        for s = 1:3
            % If range not specified, use length of time array
            ind = zeros(L,1);
            if ~exist('rnge') rnge = 1:length(ts0); end
            if ch==1
                ind(rnge) = Nind(rnge);
            elseif ch==2
                ind(rnge) = Find(rnge);
            end
            x = ts0(find(ind))';
            y = squeeze(params(p,s,find(ind)));
            omit = find( (x<=0) | isnan(x) | (y<=0) | isnan(y) );
            x(omit) = []; y(omit) = [];
            plot(x,y,'-','color',cols{s}); hold on
        end
        xlabel('time(min)','fontsize',14);
        ylabel(ylabs{p},'fontsize',14);
        title([titles{p} ' - ' chans{ch} ' ' suffix(2:end)]);
        legend('joint','mito','cyto')
        set(h(p),'color','w')
        set(gca,'fontsize',12);
        set(gcf,'PaperPositionMode','auto')
        axis tight
        saveas(h(p),[datpath 'JointMitoCytoPlots_' savelabs{p} chans{ch} '_' suffix '.fig']);
        saveas(h(p),[datpath 'JointMitoCytoPlots_' savelabs{p} chans{ch} '_' suffix '.jpg'],'jpg')
    end
    
end





