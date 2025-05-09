    function FLIM_batch_averages(path,Zrange,trange,ExcludeFolders,Label)
% Once everything has been calculated for a big batch, use this to search
% the folders for fits, take averages, then plot the average values as a
% function of sample (egg# if egg data), channel. Average over Zrange
% just have one fit file per decay, not one for NADH and one for FAD.
% -Zrange and trange: indexes of the z slices and t frames to be included
%    in the averaging. NOTE: indices start at 1, not 0, so '1' means 'first
%    z plane or time point'. Enter a cell of indices if different for different
%    samples. Enter -1 if using full range.
% -ExcludeFolders: Sometimes you don't want to analyze all folders with
%    fits in them. ExcludeFolders is a cell of folder names to exclude.
% -Label: if you want to create a custom label for the output averages file

% Version updates:
% 2018-02-20: Included provisions to handle multiple segments, 'mito, cyto,
%   and joint'. Saves to multiple pages of spreadsheet
% 2017-12-18: Incorporated 'changepoints' feature, which shows when embryos
% % undergo three developmental changes: 1st division, 2nd division (4-cell
% completion), and blastocoel opening. Specified in frames. If present, in
% will load these values into the same spreadsheet in 3 additional columns.
% 2015-05-28: Added 'BadFrs' to throw out frames that had a bad fit. Also
%  deleted all code for bead standards, which never worked. Also commented
%  out Z indexing added just previously. Uncomment if needed.
% 2015-04-21: Z indexing discrepancy. Corrected
% 2014-12-03: Adapt to include illumination profile calibration image.
% Scale intensity images to get scaled irradiance
% 2014-11-18: Update after analyzing C elegans data. Include plotting
% user-specified sub-species.
% 2014-10-27: -repair t and z plot-s
% -Divide by scan number to calc irradiances
% -add option to look at irradiance with or without std bead normalization
% Before: -constructed filenames array to order by sample, then position
% -FAD short lifetime is actually bound, so big plot changes to plot frac
% in species 1 (short lifetime)

% clear all;
% path = 'C:\Users\Tim\Documents\Academic - Research\Data\Blast_zscan\';
% % Zrange = 1:4;
% % trange = {1:2,1:4,1};

if path(end)~='\' path = [path '\']; end;
if ~exist('ExcludeFolders')|~iscell(ExcludeFolders) clear ExcludeFolders; end % ExcludeFolders must be a cell of paths, even if there is only one element
if ~exist('Label')|Label==-1 clear Label; end
close all;

sams = {};
Dfits = [];
fnames = {};
pos = {};
mask = {};
labels = {};

% Do a quick for loop. Go folder by folder (presumably sample by sample),
% get positions and sort by position. Build a 'Dfits' structure in that
% order.
Df = dir(path); Df(1:2)=[]; Df(~[Df.isdir])=[];
if exist('ExcludeFolders')
    Exc=[];
    for i = 1:length(Df)
        for j = 1:length(ExcludeFolders)
            if strcmp(Df(i).name,ExcludeFolders{j})
                Exc = [Exc i];
            end
        end
    end
    Df(Exc)=[];
end

for i = 1:length(Df)
    Dfits_tmp = dir([path Df(i).name '\*fits*.mat']);
    if isempty(Dfits_tmp)
        continue;
    end
    clear sams_tmp fnames_tmp pos_tmp mask_tmp labels_tmp
    for j = 1:length(Dfits_tmp)
        sams_tmp{j} = Df(i).name;
        fnames_tmp{j} = Dfits_tmp(j).name;
        dashes = strfind(fnames_tmp{j},'_');
        posind = strfind(fnames_tmp{j},'Pos');
        %         maskind = strfind(fnames_tmp{j},'mask');
        MasLab = fnames_tmp{j}(dashes(2)+1:strfind(fnames_tmp{j},'.mat')-1);
        
        if ~isempty(posind)
            pos_tmp{j} = (fnames_tmp{j}(posind+3:dashes(2)-1));
            %             mask_tmp{j} = fnames_tmp{j}(maskind+4:dashes(end)-1);
            mask_tmp{j} = MasLab;
            if isempty(strfind(fnames_tmp{j},'thr'))
                labels_tmp{j} = [sams_tmp{j} 'P' pos_tmp{j} mask_tmp{j}];
            else
                labels_tmp{j} = [sams_tmp{j} 'P' pos_tmp{j} 'thr'];
            end
        else
            pos_tmp{j} = -1;
            mask_tmp{j} = -1;
            labels_tmp{j} = [fnames_tmp{j}];
        end
    end
    % Sort labels
    [a,b] = sort(pos_tmp);
    Dfits_tmp = Dfits_tmp(b); pos_tmp = pos_tmp{b}; mask_tmp = mask_tmp(b); labels_tmp = labels_tmp(b);
    
    
    sams = [sams sams_tmp];
    Dfits = [Dfits; Dfits_tmp];
    fnames = [fnames fnames_tmp];
    pos = [pos pos_tmp];
    mask = [mask mask_tmp];
    labels = [labels labels_tmp];
end
usams = unique(sams);
% if no eggs, all pixels are binned with a threshold, so redox
% values can be easily calculated here. No _RedoxVals.mat file
% necessary
EggNums(1:length(Dfits)) = -1;

% For spreadsheets (one for means, one for std's)
spr{1,1} = 'Sample'; spr_std{1,1} = 'Sample';
spr{1,2} = 'NADirr'; spr_std{1,2} = 'NADirr';
spr{1,3} = 'FADirr'; spr_std{1,3} = 'FADirr';
spr{1,4} = 'RedoxRat'; spr_std{1,4} = 'RedoxRat';
spr{1,5} = 'NADtau1'; spr_std{1,5} = 'NADtau1';
spr{1,6} = 'NADtau2'; spr_std{1,6} = 'NADtau2';
spr{1,7} = 'NADfracbound'; spr_std{1,7} = 'NADfracbound';
spr{1,8} = 'FADtau1'; spr_std{1,8} = 'FADtau1';
spr{1,9} = 'FADtau2'; spr_std{1,9} = 'FADtau2';
spr{1,10} = 'FADfracbound'; spr_std{1,10} = 'FADfracbound';
spr{1,11} = 'NADtaumean'; spr_std{1,11} = 'NADtaumean';
spr{1,12} = 'FADtaumean'; spr_std{1,12} = 'FADtaumean';

% In the averages spreadsheet, use the following fields for manual labeling
spr{1,13,1} = 'Group';
spr{1,14,1} = 'GroupName';
spr{1,15,1} = 'Exclude(enter 1)';
spr{1,16,1} = 'Point Text';
spr{1,17,1} = '1st Div';
spr{1,18,1} = '2nd Div';
spr{1,19,1} = 'Blastocoel';

% Create additional dimensions for additional segments, mito, ctyo, joint
spr(:,:,2) = spr(:,:,1);
spr(:,:,3) = spr(:,:,1);


for i =1:size(Dfits,1)
    clear Ndecays Ntau1s Ntau1s_std Ntau2s Ntau2s_std Nfbound Nfbound_std Npix Nirrs Nirrs_std...
        Fdecays Ftau1s Ftau1s_std Ftau2s Ftau2s_std Ffbound Ffbound_std Fpix Firrs Firrs_std ts Zs
    disp([num2str(i) '/' num2str(size(Dfits,1))]);
    % Load fits and get param vals. Enter them into z x t sized matrice
    load([path sams{i} '\' Dfits(i).name]) % loads decays_fits_struct
    try     load([path sams{i} '\' 'multiD_indices.mat']); catch     load([path sams{i} '\' 'name_indexes.mat']); end
    
    timebins = 0;
    nonemptind = ~cellfun('isempty',decays_fits_struct); nonemptind  = find(nonemptind);
    zname = decays_fits_struct{nonemptind(1)}.filename; zdashes = strfind(zname,'_');
    %     ts0 = str2num(zname(zdashes(1)+1:zdashes(2)-1));
    %     Zs0 = str2num(zname(zdashes(3)+1:zdashes(3)+4));
    BadFrs=[];
    

    
    for seg = 1:size(decays_fits_struct{nonemptind(1)}.decay,2)% Loop over segments
        
        for j = 1:size(decays_fits_struct,1)
            if ~isempty(decays_fits_struct{j})
                restab = decays_fits_struct{j}.fit_result(:,:,seg);
                zname = decays_fits_struct{j}.filename;
                zdashes = strfind(zname,'_');
                %             ts(j) = str2num(zname(zdashes(1)+1:zdashes(2)-1))+1;
                ts(j) = nameinds{j,3}+1;
                %             if ts0==0 ts(j)=ts(j)-ts0+1; end % If t must start at 1 instead of 0
                %             Zs(j) = str2num(zname(zdashes(3)+1:zdashes(3)+4))+1;
                Zs(j) = nameinds{j,5}+1;
                %             if Zs0==0 Zs(j)=Zs(j)-Zs0+1; end % If Z must start at 1 instead of 0
                
                % Flag bad fits. Sometimes the fit doesn't work and you get extreme param
                %  values, accompanied by a 0 stddev. Check for these and flag
                %  to throw this frame out.
                BadFr = ~isempty(find(restab(:,1)>0&restab(:,2)==0&restab(:,3)==0));
                if BadFr BadFrs(Zs(j),ts(j)) = 1; end
                
                decay = decays_fits_struct{j}.decay(:,seg);
                if length(decay)>timebins timebins = length(decay); end;
                fit_start = decays_fits_struct{j}.fit_region(1); fit_end = decays_fits_struct{j}.fit_region(2);
                noise_region_from = decays_fits_struct{j}.noise_region(1);noise_region_to = decays_fits_struct{j}.noise_region(2);
                time = decays_fits_struct{j}.time;
                tstamps(Zs(j),ts(j)) = nameinds{j,8};
                % Have one param matrix for each channel (NADH vs FAD)
                if ~isempty(strfind(zname,'NADH'))
                    Ndecays{Zs(j),ts(j)} = decay;
                    Ntau1s(Zs(j),ts(j)) = restab(3,1);
                    Ntau1s_std(Zs(j),ts(j)) = restab(3,2);
                    Ntau2s(Zs(j),ts(j)) = restab(3,1)*restab(5,1);
                    Ntau2s_std(Zs(j),ts(j)) = restab(3,1)*restab(5,2);
                    Nfbound(Zs(j),ts(j)) = restab(4,1);
                    Nfbound_std(Zs(j),ts(j)) = restab(4,2);
                    
                    % Calculate irradiance values and total number of
                    % emmitters (number species 1 + 2). The latter,
                    % calculate by finding the amplitude of the exponent at
                    % t=0 (unnormalized)
                    ind = find(decays_fits_struct{j}.selected_pixel);
                    % Npix is all the pixel intensity values of the
                    % 'selected_pixels', reshaped into a linear array
                    Npix{Zs(j),ts(j)} = double(decays_fits_struct{j}.image(ind));
                    if isfield(decays_fits_struct{j},'irrSc')
                        % Check if there was a 'scaled' irr using IllProf.
                        % Otherwise, assume 'irr' was scaled. We also do this
                        Nirrs(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc(seg);
                        Nirrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc_std(seg);
                    else
                        Nirrs(Zs(j),ts(j)) = decays_fits_struct{j}.irr(seg);
                        Nirrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irr_std(seg);
                    end
                elseif ~isempty(strfind(zname,'FAD'))|~isempty(strfind(zname,'UserChan')) % If it's a user channel, just plot the irrs in the FAD channel
                    Fdecays{Zs(j),ts(j)} = decay;
                    Ftau1s(Zs(j),ts(j)) = restab(3,1);
                    Ftau1s_std(Zs(j),ts(j)) = restab(3,2);
                    Ftau2s(Zs(j),ts(j)) = restab(3,1)*restab(5,1);
                    Ftau2s_std(Zs(j),ts(j)) = restab(3,1)*restab(5,2);
                    Ffbound(Zs(j),ts(j)) = restab(4,1);
                    Ffbound_std(Zs(j),ts(j)) = restab(4,2);
                    ind = find(decays_fits_struct{j}.selected_pixel);
                    Fpix{Zs(j),ts(j)} = double(decays_fits_struct{j}.image(ind));
                    if isfield(decays_fits_struct{j},'irrSc')
                        Firrs(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc(seg);
                        Firrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc_std(seg);
                    else
                        Firrs(Zs(j),ts(j)) = decays_fits_struct{j}.irr(seg);
                        Firrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irr_std(seg);
                    end
                else
                    error('Channel error')
                end
            end
        end
        % THROW OUT bad fits
        [badZ,badt] = find(BadFrs);
        if ~isempty(badZ)
            for b = 1:length(badZ)
                Ndecays{badZ(b),badt(b)} = [];Ntau1s(badZ(b),badt(b)) = 0;Ntau1s_std(badZ(b),badt(b)) = 0;
                Ntau2s(badZ(b),badt(b)) = 0;Ntau2s_std(badZ(b),badt(b)) = 0;Nfbound(badZ(b),badt(b)) = 0;
                Nfbound_std(badZ(b),badt(b)) = 0;Npix{badZ(b),badt(b)} = [];Nirrs(badZ(b),badt(b)) = 0;Nirrs_std(badZ(b),badt(b)) = 0;
                Fdecays{badZ(b),badt(b)} = [];Ftau1s(badZ(b),badt(b)) = 0;Ftau1s_std(badZ(b),badt(b)) = 0;
                Ftau2s(badZ(b),badt(b)) = 0;Ftau2s_std(badZ(b),badt(b)) = 0;Ffbound(badZ(b),badt(b)) = 0;
                Ffbound_std(badZ(b),badt(b)) = 0;Fpix{badZ(b),badt(b)} = [];Firrs(badZ(b),badt(b)) = 0;Firrs_std(badZ(b),badt(b)) = 0;
            end
        end
        
        % Calculate dt (in minutes) for plotting later
        % Find z planes that found masks for the first and last time point
        % First get rid of frames that were entirely black
        tstamps2 = tstamps; RemCl = [];
        for cl = 1:size(tstamps2,2)
            if isempty(find(tstamps2(:,cl))) RemCl = [RemCl cl]; end
        end
        tstamps2(:,RemCl)=[];
        ind = intersect(find(tstamps2(:,1)),find(tstamps2(:,end)));
        tmp = tstamps2(ind(1),:); tmp = tmp(find(tmp));
        dts{i} = (tmp(end)-tmp(1))/length(tmp)*86400/60;
        
        % Reshape params to z x t matrices. First pad with zeros
        SmInd = min([size(nameinds,1),size(decays_fits_struct,1)]);
        nametmp = nameinds(1:SmInd,:); decaystmp = decays_fits_struct(1:SmInd,:);
        tgrid = unique(sort([nametmp{( ([nametmp{:,7}]~=-1) & (~cellfun('isempty',decaystmp))' ),3}]));
        tgrid = tgrid + 1; %tgrid - min(tgrid) + 1; % Not sure why I was subtracting the minimum
        Zgrid = unique(sort([nametmp{[nametmp{:,7}]~=-1,5}]));
        Zgrid = Zgrid + 1; % - min(Zgrid)
        % If range cells entered:
        rangeind = find(strcmp(usams,sams{i}));
        if exist('Zrange')
            if iscell(Zrange)
                if Zrange{rangeind}~=-1
                    samZrange = Zrange{rangeind};
                else
                    samZrange = Zgrid;
                end
            elseif Zrange~=-1
                samZrange = Zrange;
            else
                samZrange = Zgrid;
            end
        else
            samZrange = Zgrid;
        end
        if exist('trange')
            if iscell(trange)
                if trange{rangeind}~=-1
                    samtrange = trange{rangeind};
                else
                    samtrange = tgrid;
                end
            elseif trange~=-1
                samtrange = trange;
            else
                samtrange = tgrid;
            end
        else
            samtrange = tgrid;
        end
        if length(samZrange)>length(Zgrid) samZrange = Zgrid; end;
        if length(samtrange)>length(tgrid) samtrange = tgrid; end;
        % If any of the time points for this
        
        % Make tau1's the short lifetime
        if exist('Ntau1s')
            if size(Ntau1s,2)<length(tgrid)
                % If not enough rows/columns, pad with 0's to get the same
                % dimensions as Zgrid and tgrid
                Ntau1s(:,end+1:length(tgrid))=0; Ntau1s_std(:,end+1:length(tgrid))=0;
                Ntau2s(:,end+1:length(tgrid))=0; Ntau2s_std(:,end+1:length(tgrid))=0;
                Nfbound(:,end+1:length(tgrid))=0; Nfbound_std(:,end+1:length(tgrid))=0;
                Npix(:,end+1:length(tgrid))=num2cell(0);
                Nirrs(:,end+1:length(tgrid))=0; Nirrs_std(:,end+1:length(tgrid))=0;
                
            end
            if size(Ntau1s,1)<length(Zgrid)
                Ntau1s(end+1:length(Zgrid),:)=0; Ntau1s_std(end+1:length(Zgrid),:)=0;
                Ntau2s(end+1:length(Zgrid),:)=0; Ntau2s_std(end+1:length(Zgrid),:)=0;
                Nfbound(end+1:length(Zgrid),:)=0; Nfbound_std(end+1:length(Zgrid),:)=0;
                Npix(end+1:length(Zgrid),:)=num2cell(0);
                Nirrs(end+1:length(Zgrid),:)=0; Nirrs_std(end+1:length(Zgrid),:)=0;
            end
        end
        if exist('Ftau1s')
            if size(Ftau1s,2)<length(tgrid)
                Ftau1s(:,end+1:length(tgrid))=0; Ftau1s_std(:,end+1:length(tgrid))=0;
                Ftau2s(:,end+1:length(tgrid))=0; Ftau2s_std(:,end+1:length(tgrid))=0;
                Ffbound(:,end+1:length(tgrid))=0; Ffbound_std(:,end+1:length(tgrid))=0;
                Fpix(:,end+1:length(tgrid))=num2cell(0);
                Firrs(:,end+1:length(tgrid))=0; Firrs_std(:,end+1:length(tgrid))=0;
            end
            if size(Ftau1s,1)<length(Zgrid)
                Ftau1s(end+1:length(Zgrid),:)=0; Ftau1s_std(end+1:length(Zgrid),:)=0;
                Ftau2s(end+1:length(Zgrid),:)=0; Ftau2s_std(end+1:length(Zgrid),:)=0;
                Ffbound(end+1:length(Zgrid),:)=0; Ffbound_std(end+1:length(Zgrid),:)=0;
                Fpix(end+1:length(Zgrid),:)=num2cell(0);
                Firrs(end+1:length(Zgrid),:)=0; Firrs_std(end+1:length(Zgrid),:)=0;
            end
        end
        
        
        % AVERAGES
        % AVERAGES CALCULATED WITH WEIGHTING FACTORS = 1/SIGMA^2
        % http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
        
        % Average over Z scans, keep time points
        % Use only non-zero values, where params were successfully calced
        for k = samtrange % Average one column at a time to exclude zeros
            if exist('Ntau1s')
                tmp = Ntau1s(samZrange,k); tmp_std = Ntau1s_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ntau1_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ntau1_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Ntau2s(samZrange,k); tmp_std = Ntau2s_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ntau2_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ntau2_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Nfbound(samZrange,k); tmp_std = Nfbound_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Nfbound_Z_aves{i}(k) =  1 - sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Nfbound_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Nirrs(samZrange,k); tmp_std = Nirrs_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Nirr_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Nirr_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
            end
            
            if exist('Ftau1s')
                tmp = Ftau1s(samZrange,k); tmp_std = Ftau1s_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ftau1_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ftau1_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Ftau2s(samZrange,k); tmp_std = Ftau2s_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ftau2_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ftau2_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Ffbound(samZrange,k); tmp_std = Ffbound_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ffbound_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ffbound_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Firrs(samZrange,k); tmp_std = Firrs_std(samZrange,k); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Firr_Z_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Firr_Z_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
            end
        end
        
        if exist('Ntau1s')&exist('Ftau1s')
            redoxirr_Z_aves{i} = Nirr_Z_aves{i}./Firr_Z_aves{i};
            redoxirr_Z_aves_std{i} = Nirr_Z_aves{i}./Firr_Z_aves{i}.*sqrt((Nirr_Z_aves_std{i}./Nirr_Z_aves{i}).^2+(Firr_Z_aves_std{i}./Firr_Z_aves{i}).^2);
        end
        % Average over t points, keep z resolution and see how params vary
        % with z
        for k = samZrange % Average one column at a time to exclude zeros
            if exist('Ntau1s')
                tmp = Ntau1s(k,samtrange); tmp_std = Ntau1s_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ntau1_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ntau1_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Ntau2s(k,samtrange); tmp_std = Ntau2s_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ntau2_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ntau2_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Nfbound(k,samtrange); tmp_std = Nfbound_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Nfbound_t_aves{i}(k) =  1 - sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Nfbound_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Nirrs(k,samtrange); tmp_std = Nirrs_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Nirr_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Nirr_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
            end
            
            if exist('Ftau1s')
                tmp = Ftau1s(k,samtrange); tmp_std = Ftau1s_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ftau1_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ftau1_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Ftau2s(k,samtrange); tmp_std = Ftau2s_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ftau2_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ftau2_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Ffbound(k,samtrange); tmp_std = Ffbound_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Ffbound_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Ffbound_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
                tmp = Firrs(k,samtrange); tmp_std = Firrs_std(k,samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
                Firr_t_aves{i}(k) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
                Firr_t_aves_std{i}(k) = sqrt(1./sum(1./tmp_std.^2));
            end
        end
        
        if exist('Ntau1s')&exist('Ftau1s')
            redoxirr_t_aves{i} = Nirr_t_aves{i}./Firr_t_aves{i};
            redoxirr_t_aves_std{i} = Nirr_t_aves{i}./Firr_t_aves{i}.*sqrt((Nirr_t_aves_std{i}./Nirr_t_aves{i}).^2+(Firr_t_aves_std{i}./Firr_t_aves{i}).^2);
        end
        
        % Sometimes short and long lifetimes get switched, but I want tau1
        % to be the SHORT lifetime, so switch if this happens
        if exist('Ntau1s')
            if mean(Ntau2_t_aves{i}(~isnan(Ntau2_t_aves{i})))<mean(Ntau1_t_aves{i}(~isnan(Ntau1_t_aves{i})))
                tmp = Ntau2_Z_aves{i}; Ntau2_Z_aves{i} = Ntau1_Z_aves{i}; Ntau1_Z_aves{i} = tmp;
                tmp = Ntau2_Z_aves_std{i}; Ntau2_Z_aves_std{i} = Ntau1_Z_aves_std{i}; Ntau1_Z_aves_std{i} = tmp;
                tmp = Ntau2_t_aves{i}; Ntau2_t_aves{i} = Ntau1_t_aves{i}; Ntau1_t_aves{i} = tmp;
                tmp = Ntau2_t_aves_std{i}; Ntau2_t_aves_std{i} = Ntau1_t_aves_std{i}; Ntau1_t_aves_std{i} = tmp;
                Nfbound_Z_aves{i} = 1 - Nfbound_Z_aves{i}; Nfbound_t_aves{i} = 1 - Nfbound_t_aves{i};
            end
        end
        if exist('Ftau1s')
            if mean(Ftau2_t_aves{i}(~isnan(Ftau2_t_aves{i})))<mean(Ftau1_t_aves{i}(~isnan(Ftau1_t_aves{i})))
                tmp = Ftau2_Z_aves{i}; Ftau2_Z_aves{i} = Ftau1_Z_aves{i}; Ftau1_Z_aves{i} = tmp;
                tmp = Ftau2_Z_aves_std{i}; Ftau2_Z_aves_std{i} = Ftau1_Z_aves_std{i}; Ftau1_Z_aves_std{i} = tmp;
                tmp = Ftau2_t_aves{i}; Ftau2_t_aves{i} = Ftau1_t_aves{i}; Ftau1_t_aves{i} = tmp;
                tmp = Ftau2_t_aves_std{i}; Ftau2_t_aves_std{i} = Ftau1_t_aves_std{i}; Ftau1_t_aves_std{i} = tmp;
                Ffbound_Z_aves{i} = 1- Ffbound_Z_aves{i}; Ffbound_t_aves{i} = 1- Ffbound_t_aves{i};
            end
        end
        
        % Average t (selected) points for a big final plot of all samples. Use
        % weights again.
        if exist('Ntau1s')
            tmp = Ntau1_Z_aves{i}(samtrange); tmp_std = Ntau1_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Ntau1aves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Ntau1aves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            tmp = Ntau2_Z_aves{i}(samtrange); tmp_std = Ntau2_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Ntau2aves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Ntau2aves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            tmp = Nfbound_Z_aves{i}(samtrange); tmp_std = Nfbound_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Nfshortaves_wt(i) =  1 - sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Nfshortaves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            tmp = Nirr_Z_aves{i}(samtrange); tmp_std = Nirr_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Nirraves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Nirraves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            % Write results in a spreadsheet to export to xls
            spr{i+1,2,seg} = Nirraves_wt(i);
            spr{i+1,5,seg} = Ntau1aves_wt(i);
            spr{i+1,6,seg} = Ntau2aves_wt(i);
            spr{i+1,7,seg} = 1-Nfshortaves_wt(i); % NADH fbound
            spr{i+1,11,seg} = Ntau1aves_wt(i)*(1-Nfshortaves_wt(i))+Ntau2aves_wt(i)*Nfshortaves_wt(i);
            % Std dev spreadsheet
            spr_std{i+1,2,seg} = Nirraves_wt_std(i);
            spr_std{i+1,5,seg} = Ntau1aves_wt_std(i);
            spr_std{i+1,6,seg} = Ntau2aves_wt_std(i);
            spr_std{i+1,7,seg} = 1-Nfshortaves_wt_std(i); % NADH fbound
            spr_std{i+1,11,seg} = Ntau1aves_wt_std(i)*(1-Nfshortaves_wt_std(i))+Ntau2aves_wt_std(i)*Nfshortaves_wt_std(i);
            
        end
        
        if exist('Ftau1s')
            tmp = Ftau1_Z_aves{i}(samtrange); tmp_std = Ftau1_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Ftau1aves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Ftau1aves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            tmp = Ftau2_Z_aves{i}(samtrange); tmp_std = Ftau2_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Ftau2aves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Ftau2aves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            tmp = Ffbound_Z_aves{i}(samtrange); tmp_std = Ffbound_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Ffshortaves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Ffshortaves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            tmp = Firr_Z_aves{i}(samtrange); tmp_std = Firr_Z_aves_std{i}(samtrange); tmp=tmp(tmp>0&tmp<Inf);tmp_std=tmp_std(tmp_std>0&tmp_std<Inf);
            Firraves_wt(i) = sum(tmp./tmp_std.^2)./sum(1./tmp_std.^2);
            Firraves_wt_std(i) = sqrt(1./sum(1./tmp_std.^2));
            % Spreadsheet
            spr{i+1,1,seg} = labels{i};
            spr{i+1,3,seg} = Firraves_wt(i);
            spr{i+1,8,seg} = Ftau1aves_wt(i);
            spr{i+1,9,seg} = Ftau2aves_wt(i);
            spr{i+1,10,seg} = Ffshortaves_wt(i); %
            spr{i+1,12,seg} = Ftau1aves_wt(i)*(1-Ffshortaves_wt(i))+Ftau2aves_wt(i)*Ffshortaves_wt(i);
            % Std dev spreadsheet
            spr_std{i+1,1,seg} = labels{i};
            spr_std{i+1,3,seg} = Firraves_wt(i);
            spr_std{i+1,8,seg} = Ftau1aves_wt(i);
            spr_std{i+1,9,seg} = Ftau2aves_wt(i);
            spr_std{i+1,10,seg} = Ffshortaves_wt(i); %
            spr_std{i+1,12,seg} = Ftau1aves_wt(i)*(1-Ffshortaves_wt(i))+Ftau2aves_wt(i)*Ffshortaves_wt(i);
        end
        
        if exist('Ntau1s')&exist('Ftau1s')
            redoxirraves_wt(i) = Nirraves_wt(i)./Firraves_wt(i);
            redoxirraves_wt_std(i) = Nirraves_wt(i)./Firraves_wt(i).*sqrt((Nirraves_wt_std(i)./Nirraves_wt(i)).^2+(Firraves_wt_std(i)./Firraves_wt(i)).^2);
            spr{i+1,4,seg} = redoxirraves_wt(i);
            spr_std{i+1,4,seg} = redoxirraves_wt_std(i);
        end
        
        % Get embryo label for plotting
        dashes = strfind(Dfits(i).name,'_');
        PosInd = strfind(Dfits(i).name,'Pos');
        MaskInd = strfind(Dfits(i).name,'mask');
        PosNum = Dfits(i).name(PosInd+3:dashes(2)-1);
        if length(dashes)>2 % E.g. '..._fxshft.mat'
            MaskNum = Dfits(i).name(MaskInd+4:dashes(3)-1);
        else
            MaskNum = Dfits(i).name(MaskInd+4:end-4);
        end
        
        % Changepoints: search and load into spr, if present
        Dcp = dir([Dfits(i).folder '\embryo_data_changepoints.mat']);
        if ~isempty(Dcp)
            load([Dcp.folder '\' Dcp.name])
            Pcpind = strfind(chngpnts_accum(:,1),['Pos' PosNum]);
            Pcpind = ~cellfun('isempty',Pcpind);
            Mcpind = strfind(chngpnts_accum(:,1),['mask' MaskNum]);
            Mcpind = ~cellfun('isempty',Mcpind);
            cpind = find(Pcpind&Mcpind);
            cps = chngpnts_accum{cpind,2};
            if length(cpind)>1
                error('Unique labeling of changepoints violated.')
            end
            % Note: the frames in chngpnts_accum are the time points, and
            % that's actually correct. FLIM_time_plots.m uses those same units
            % (ie, 'time point', not 'frame num'
            spr{i+1,17,seg} = cps(1);
            spr{i+1,18,seg} = cps(2);
            spr{i+1,19,seg} = cps(3);
        end
        
        dashes = strfind(sams{i},'_');
        if strcmp(sams{i}(1),'s')
            if isempty(dashes)
                spr{i+1,16} = [sams{i} 'P' PosNum 'm' MaskNum];
            else
                spr{i+1,16} = [sams{i}(1:dashes(1)-1) 'P' PosNum 'm' MaskNum];
            end
        else
            spr{i+1,16} = ['P' PosNum 'm' MaskNum];
        end
        
    end % Seg loop end
    
    
end

% Write averages from different segments into different worksheets of the
% excell spreadsheet. (Use reverse order so it's on the first page when you
% open up the spreadsheet)
if exist('Label') delete([path '\' Label '.xls']); delete([path '\' Label '_std.xls']); end
delete([path '\ParamsAllSams.xls']); delete([path '\ParamsAllSams_std.xls']);
for seg = size(decays_fits_struct{nonemptind(1)}.decay,2):-1:1
    if exist('Label')
        xlswrite([path '\' Label '.xls'],spr(:,:,seg),seg);
        xlswrite([path '\' Label '_std.xls'],spr_std(:,:,seg),seg);
    else
        xlswrite([path '\ParamsAllSams.xls'],spr(:,:,seg),seg);
        xlswrite([path '\ParamsAllSams_std.xls'],spr_std(:,:,seg),seg);
    end
end





