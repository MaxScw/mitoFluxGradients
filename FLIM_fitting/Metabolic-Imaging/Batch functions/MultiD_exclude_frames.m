function MultiD_exclude_frames(acqpath,exframes,poss)
% After tiff conversion, look at the tiff image sequences. If any of them
% are bad for some reason, they can be excluded before further processing.
% This function:
% -deletes the tiff
% -deletes the corresponding sdt's
% -adjusts 'nameinds'

% clear all;
% acqpath = 'C:\Users\Tim\Documents\Academic - Research\Data\REWRITE2EXP\s1_a1_Day1_test';
% exframes = 1:3;
% poss = 0;

if acqpath(end)~='\'; acqpath = [acqpath '\']; end;

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
try     load([acqpath 'multiD_indices.mat']); catch     load([acqpath 'name_indexes.mat']); end
% Load both JointMasks and default mdpars.
if exist([acqpath 'multiD_pars_JointMasks.mat']) 
    load([acqpath 'multiD_pars_JointMasks.mat']); 
    mdparsJ = mdpars;
end
if exist([acqpath 'multiD_pars.mat']) load([acqpath 'multiD_pars.mat']); end

for posnum = 1:length(Dpos)
    
    D = dir([sdtpath Dpos(posnum).name '*.sdt']);
    PosNum = Dpos(posnum).name(4:end);
    PosInd = strcmp(nameinds(:,3),PosNum); PosInd = PosInd';
    
    % Maybe you only want to do certain positions, like if you want to redo
    % certain positions with different image processing parameters
    if exist('poss')&poss~=-1
        if ~strcmp(num2str(poss),PosNum)
            continue;
        end
    end
    
    tifpath = [sdtpath 'IntTiffs_' Dpos(posnum).name '_' Run '\'];
    ROIpath = [sdtpath 'ROIsCheck_' Dpos(posnum).name '_' Run '\'];
    ROIGenpath = [sdtpath 'ROIsCheckGen_' Dpos(posnum).name '_' Run '\'];
    FLIMpath = [sdtpath 'FLIMageTiffs_' Dpos(posnum).name '_' Run '\'];
    SHGpath = [sdtpath 'SHGTiffs_' Dpos(posnum).name '_' Run '\'];
    
    for i = exframes
        if i>max([nameinds{:,6}])
            continue;
        end
        % Delete ROICheck png's
        Dt = dir([ROIpath 'ROI_im' num2str(i,'%03i') '.png']);
        if ~isempty(Dt) delete([ROIpath Dt(1).name]); end
        DtGen = dir([ROIGenpath 'ROI_im' num2str(i,'%03i') '.png']);
        if ~isempty(DtGen) delete([ROIGenpath DtGen(1).name]); end
        % Don't delete the raw data, so you can see it if you want to.
        % Check ROI movies to see what frames were deleted
        %         % Delete tiff:
        %         Dt = dir([tifpath 'fr' num2str(i,'%05i') '.tif']);
        %         if ~isempty(Dt) delete([tifpath Dt(1).name]); end
        %         Dt = dir([FLIMpath 'fr' num2str(i,'%05i') '.tif']);
        %         if ~isempty(Dt) delete([FLIMpath Dt(1).name]); end
        %         Dt = dir([SHGpath 'fr' num2str(i,'%05i') '.tif']);
        %         if ~isempty(Dt) delete([SHGpath Dt(1).name]); end
        
        % Delete sdt files
        ind = (PosInd)&([nameinds{:,6}]==i);
        if isempty(find(ind))
            continue;
        end
        t(i) = nameinds{find(ind),2};
        z(i) = nameinds{find(ind),5};
        nameinds(find(ind),[7 8]) = num2cell(-1);
        %         Dt = dir([sdtpath Dpos(posnum).name '\sdt_' num2str(t(i),'%09g') '_NADH_' num2str(z(i),'%03g') '.sdt']);
        %         if ~isempty(Dt) delete([sdtpath Dpos(posnum).name '\' Dt(1).name]); end
        %         Dt = dir([sdtpath Dpos(posnum).name '\sdt_' num2str(t(i),'%09g') '_FAD_' num2str(z(i),'%03g') '.sdt']);
        %         if ~isempty(Dt) delete([sdtpath Dpos(posnum).name '\' Dt(1).name]); end
        %         Dt = dir([sdtpath Dpos(posnum).name '\sdt_' num2str(t(i),'%09g') '_UserChan_' num2str(z(i),'%03g') '.sdt']);
        %         if ~isempty(Dt) delete([sdtpath Dpos(posnum).name '\' Dt(1).name]); end
     
        % Set elements of mdpars to NaN
        % For ref, mdpars was defined as
        % mdpars(:,:,Tind,Pind,Chind,Zind,Mind,seg)=[fit_result(:,1:2);[irr irr_stderr];[decay_struct{j}.timestp decay_struct{j}.timestp]];
        if exist([acqpath 'multiD_pars_JointMasks.mat']) mdparsJ(:,:,t(i)+1,posnum,:,z(i)+1,:,:)=NaN; end
        if exist([acqpath 'multiD_pars.mat']) mdpars(:,:,t(i)+1,posnum,:,z(i)+1,:,:)=NaN; end
    end
    % If decays and fits have been calculated, clear these entries so they
    % don't get used in the averages or plots. Compile 'mask' files
    % (individual eggs) and 'GenMask' files (single mask) into one dir
    % array.
    Dd = dir([acqpath 'decays_' Dpos(posnum).name '_mask*.mat']);
    Dd = [Dd dir([acqpath 'decays_' Dpos(posnum).name '_SingleMasks*.mat'])];
    Df = dir([acqpath 'fits_' Dpos(posnum).name '_mask*.mat']);
    Df = [Df dir([acqpath 'fits_' Dpos(posnum).name '_SingleMasks*.mat'])];
    if ~isempty(Dd)
        for fnum = 1:size(Dd,1)
            load([acqpath Dd(fnum).name]);
            if ~isempty(Df) load([acqpath Df(fnum).name]); end
            for i = exframes
                if i<=max([nameinds{:,6}])
                    decaysind = [nameinds{:,3}]==PosNum&[nameinds{:,2}]==t(i)&[nameinds{:,5}]==z(i);
                    decay_struct(decaysind)= cell(1);
                    if ~isempty(Df) decays_fits_struct(decaysind)= cell(1); end
                end
            end
            save([acqpath Dd(fnum).name],'decay_struct')
            if ~isempty(Df) save([acqpath Df(fnum).name],'decays_fits_struct'); end
        end
    end
    
    % If the whole position was removed, remove the folders as well
    ind = (PosInd)&([nameinds{:,7}]>-1);
    if isempty(find(ind))
        [a,b,c] = rmdir([sdtpath Dpos(posnum).name],'s');
        [a,b,c] = rmdir(tifpath,'s');
        [a,b,c] = rmdir(ROIpath,'s');
        [a,b,c] = rmdir(FLIMpath,'s');
        [a,b,c] = rmdir(SHGpath,'s');
    end
    
end
save([acqpath 'multiD_indices.mat'],'nameinds');
if exist([acqpath 'multiD_pars.mat']) save([acqpath 'multiD_pars.mat'],'mdpars'); end
if exist([acqpath 'multiD_pars_JointMasks.mat']) 
    mdpars = mdparsJ;
    save([acqpath 'multiD_pars_JointMasks.mat'],'mdpars');
end