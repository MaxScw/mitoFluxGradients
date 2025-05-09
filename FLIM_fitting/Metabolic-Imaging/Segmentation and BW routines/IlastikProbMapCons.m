function IlastikProbMapCons(acqpath,SaveBool)
% This function reads in the probability maps created by
% IlastikProbMap_shell.m and applies a dynamic adjustment to them to
% conserve the total probabily density of each segment (mito, cyto). It
% matches each map to the average of the first 5 time points.
% INPUTS:
% -acqpath: path to acquisition data
% -SaveBool: '1' means save the adjusted prob maps in an additional folder,
%   'ProbMaps... adjusted'
% OUTPUTS: Saves a simple array, 'ProbMapExps', of exponents that best
% match the intended prob density. This list of exponents is used later by
% FLIM_decay_from_probmap.m

% clear all; close all;
% acqpath = 'Z:\Lab\Emily\2016-11-10 Emily O2 Drop on Embryos\2017-06-01 1-cell NO O2 DROP\s1_a1_1cell\';
% SaveBool = 1;

% Data paths
if acqpath(end)~='\'; acqpath = [acqpath '\']; end;
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end
NADHBool = 1; FADBool = 1;
if isempty(cell2mat(strfind(nameinds(:,4),'NADH'))) NADHBool = 0; end
if isempty(cell2mat(strfind(nameinds(:,4),'FAD'))) FADBool = 0; end
if ~exist('SaveBool')|SaveBool==-1 SaveBool = 0; end
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

if ~isempty(Dpos)
    srtpath = [acqpath 'sorted_sdts\'];
    for i = 1:size(Dpos,1)
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
        Probpaths{i} = [srtpath 'ProbMaps_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir([Probpaths{i}]);
        if SaveBool [a,b] = mkdir([Probpaths{i}(1:end-1) '_adjusted']); end
        
        %         % save figures with overlaid ROIs to do quick checks after batch processing
        %         ROIpaths{i} = [srtpath 'ROIsCheckGen_' Dpos(i).name '_' Run '\'];
        %         [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = path;
end

for posnum = 1:size(Dpos,1)
    % posnum =1;
    uManPos = Dpos(posnum).name(4:end);
    %     strnums = sscanf(uManPos ,'%g'); %Find the numbers in the name
    %     uManPos = strnums(1); % Assume name starts with 'Pos#' and the first number is the pos number
    %     PosInd = strcmp(nameinds(:,2),uManPos); PosInd = PosInd';
    
    %     % Maybe you only want to do certain positions, like if you want to redo
    %     % certain positions with different image processing parameters
    %     if exist('poss')&poss~=-1
    %         if ~strcmp(num2str(poss),uManPos)
    %             continue;
    %         end
    %     end
    
    %     frames = unique(sort([nameinds{PosInd&[nameinds{:,7}]>-1,6}]));
    %     frames = frames(frames>0);
    %     ts = unique(sort([nameinds{PosInd,3}]))+1;
    %     Zs = unique(sort([nameinds{PosInd,5}]))+1;
    
    % +1 because uMan indexes start at 0, I like 1
    if exist('Frs')&Frs~=-1 frames = Frs; end;
    
    set(gcf,'PaperPositionMode', 'auto')
    clear Nmasks Fmasks Masks EggVals NADim FADim
    
    % Only code for the NADH-FAD product segmentation. This is what we
    % typically use for perturbation experiments.
    Dprb = dir([Probpaths{posnum} '\*_P*.tif']);
    im1=imread([Probpaths{1} Dprb(1).name]);
    [ydim,xdim] = size(im1);
    Probs = zeros([size(im1) length(Dprb)]);
    
    % First, get total prob densities for all channel and frames
    for fr = 1:length(Dprb)
        Prob=imread([Probpaths{posnum} Dprb(fr).name]);
        
        Tots = squeeze(sum(sum(Prob,2),1));
        % Find background
        BG = Prob(:,:,Tots==max(Tots));
        % Get index for forground segments, to skip calculations for
        % background
        %         FGind = 1:size(Prob,3); FGind(Tots==max(Tots))=[];
        for seg = 1:size(Prob,3)
            segprob = Prob(:,:,seg);
            % Only include probability mass in area that is not BG
            if seg~=find(Tots==max(Tots)) segprob(BG>100) = 0; end
            SegProbTots(seg,fr) = sum(sum(segprob));
            Probs(:,:,seg,fr) = segprob;
        end
        
    end
    
    % Get average of first 5 frames
    SegProbTots5ave = mean(SegProbTots(:,1:5),2);
    
    % Plot to see if prob densities are changing over experiment
    h=figure;
    frs = 1:size(SegProbTots,2);
    plot(frs,SegProbTots(1,:)./SegProbTots5ave(1),'r',frs,SegProbTots(2,:)./SegProbTots5ave(2),'g'...
        ,frs,SegProbTots(3,:)./SegProbTots5ave(3),'b',frs,SegProbTots(4,:)./SegProbTots5ave(4),'k')
    legend('mito','cyto','bg','nuc');
    xlabel('frame'); ylabel('Rel change in prob mass')
    saveas(h,[Probpaths{posnum}(1:end-1) '_adjusted\TotalProbOverExp.fig'])
    close(h)
    
    % Now go through and adjust all the prob maps by exponentiating them.
    for fr = 1:length(Dprb)
        disp(['Frame: ' num2str(fr)]);
        % look at a range of exponent values and see what gives the closest
        % prob mass to the 5-fr ave
        erng = .5:.01:2;
        clear BestProb
        
        for seg = 1:size(Prob,3)
            for e = 1:length(erng)
                ExpProb = Probs(:,:,seg,fr).^erng(e); ExpProb = ExpProb./max(ExpProb(:)).*255;
                ExpTots(e) = sum(sum(ExpProb ));
%                 % Optional imwrite
%                 if fr==1
%                     imwrite(uint8(ExpProb),[srtpath '\fr1_exp_rng\fr1_' num2str(erng(e)*100,'%3g') '.tif'])
%                 end
            end
            TotDif = abs(ExpTots-SegProbTots5ave(seg));
            % Optional plot to show how clear the
            %         semilogy(erng,TotDif); xlabel('Exponent'); ylabel('Diff w/ 5-fr ave');
            BestInd = find(TotDif==min(TotDif));
            
            % Load into the frame that corresponds with the frame number in
            % the tiff's file name. 
            AcqFr = str2num(Dprb(fr).name(3:7));
            BestExps(seg,AcqFr) = erng(BestInd(1));
            ExpProb = Probs(:,:,seg,fr).^BestExps(seg,fr); ExpProb = ExpProb./max(ExpProb(:)).*255;
            BestProb(:,:,seg) = ExpProb;
            AdjustedTots(seg,fr) = sum(sum(ExpProb));
        end
        
        if SaveBool
            imwrite(uint8(BestProb(:,:,1:3)),[Probpaths{posnum}(1:end-1) '_adjusted\' Dprb(fr).name])
        end
        
    end
    save([acqpath 'ProbBestExps_' Dpos(posnum).name '.mat'],'BestExps')
end


% Plot of best exponents for all frames
h=figure;
frs = 1:size(BestExps,2);
plot(frs,BestExps(1,:),'r',frs,BestExps(2,:),'g'...
    ,frs,BestExps(3,:),'b',frs,BestExps(4,:),'k')
legend('mito','cyto','bg','nuc');
xlabel('frame'); ylabel('Optimal exponent for matching prob mass')
saveas(h,[Probpaths{posnum}(1:end-1) '_adjusted\BestExps.fig'])
close(h)

% % Plot adjusted totals (temp sanity check)
% figure; 
% h=figure;
% frs = 1:size(SegProbTots,2);
% plot(frs,AdjustedTots(1,:),'r',frs,AdjustedTots(2,:),'g'...
%     ,frs,AdjustedTots(3,:),'b',frs,AdjustedTots(4,:),'k')
% legend('mito','cyto','bg','nuc');
% xlabel('frame'); ylabel('Prob mass of adjusted frames')






