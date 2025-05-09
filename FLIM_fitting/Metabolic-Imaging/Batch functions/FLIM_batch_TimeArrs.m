function FLIM_batch_TimeArrs(path,Zrange,trange,ExcludeFolders,Label,tstmpBool)
% Once everything has been calculated for a big batch, use this to search
% the folders for fits, take averages, then plot the average values as a
% function of sample (egg# if egg data), channel. Average over Zrange
% just have one fit file per decay, not one for NADH and one for FAD.
% INPUTS:
% -Zrange and trange: indexes of the z slices and t frames to be included
%    in the averaging. Enter a cell of indices if different for different
%    samples. Enter -1 if using full range.
% -ExcludeFolders: Sometimes you don't want to analyze all folders with
%    fits in them. ExcludeFolders is a cell of folder names to exclude.
% -tstmpBool: Enter 1 if you want the time arrays to retain the absolute
%    timestamp values of each time point. Default (0) is to call the first
%    time point or each mask trajectory 't=0'
% -Label: if you want to create a custom label for the output TimeArrs file

% Version updates:
% 2018-02-20: Included provisions to handle multiple segments, 'mito, cyto,
%   and joint'. Saves to same TimeArrs.mat file, but TimeArrs struct has an
%   additional dimension to it for the multiple segments
% 2017-09-14: Cleaned up code by removing ..._t_aves variables, since they
%  aren't used. Also added tstmpBool option to retain absolute time stamps.
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
% path = 'C:\Users\Tim\Documents\Academic - Research\Data\2017-05-19 Batch 2\';
% % ExcludeFolders = {'2016-04-06 Repeat','FAD_int_DualCh','FAD_int_tests_cross_750nm','FAD_int_tests_cross_845nm'};
% Zrange = -1;
% trange = -1;

if path(end)~='\' path = [path '\']; end;
if ~exist('tstmpBool')|tstmpBool==-1 tstmpBool = 0; end
if ~exist('ExcludeFolders')|~iscell(ExcludeFolders) clear ExcludeFolders; end % ExcludeFolders must be a cell of paths, even if there is only one element
if ~exist('Label')|Label==-1 clear Label; end
close all;

sams = {};
Dfits = [];
fnames = {};
pos = [];
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
                labels_tmp{j} = [sams_tmp{j} 'P' num2str(pos_tmp{j}) mask_tmp{j}];
            else
                labels_tmp{j} = [sams_tmp{j} 'P' num2str(pos_tmp{j}) 'thr'];
            end
        else
            pos_tmp{j} = -1;
            mask_tmp{j} = -1;
            labels_tmp{j} = [fnames_tmp{j}];
        end
    end
    % Sort labels
    [a,b] = sort(pos_tmp);
    Dfits_tmp = Dfits_tmp(b); pos_tmp = pos_tmp(b); mask_tmp = mask_tmp(b); labels_tmp = labels_tmp(b);
    
    
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



for i =1:size(Dfits,1)
    clear Ndecays Ntau1s Ntau1s_std Ntau2s Ntau2s_std Nfbound Nfbound_std Npix Nirrs Nirrs_std...
        Fdecays Ftau1s Ftau1s_std Ftau2s Ftau2s_std Ffbound Ffbound_std Fpix Firrs Firrs_std ts Zs...
        Ntimestp Ftimestp
    disp([num2str(i) '/' num2str(size(Dfits,1))]);
    % Load fits and get param vals. Enter them into z x t sized matrice
    load([path sams{i} '\' Dfits(i).name]) % loads decays_fits_struct
    try     load([path sams{i} '\' 'multiD_indices.mat']); catch     load([path sams{i} '\' 'name_indexes.mat']); end
    
    % Load irf's for each sample
    irfind = find(~cellfun('isempty',decays_fits_struct));
    irfs{i} = decays_fits_struct{irfind(1)}.irf;
    
    timebins = 0;
    nonemptind = ~cellfun('isempty',decays_fits_struct); nonemptind  = find(nonemptind);
    zname = decays_fits_struct{nonemptind(1)}.filename; zdashes = strfind(zname,'_');
    %     ts0 = str2num(zname(zdashes(1)+1:zdashes(2)-1));
    %     Zs0 = str2num(zname(zdashes(3)+1:zdashes(3)+4));
    BadFrs=[];
    
    for seg = 1:size(decays_fits_struct{nonemptind(1)}.decay,2) % Loop over segments
        
        for j = 1:size(decays_fits_struct,1)
            if ~isempty(decays_fits_struct{j})
                restab = decays_fits_struct{j}.fit_result(:,:,seg);
                zname = decays_fits_struct{j}.filename;
                zdashes = strfind(zname,'_');
                
                ts(j) = nameinds{j,3}+1;
                Zs(j) = nameinds{j,5}+1;
                
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
                if strfind(zname,'NADH')
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
                    Nnumpix(Zs(j),ts(j)) = length(Npix{Zs(j),ts(j)});
                    Nphot(Zs(j),ts(j)) = sum(Npix{Zs(j),ts(j)});
                    if isfield(decays_fits_struct{j},'irrSc')
                        % Check if there was a 'scaled' irr using IllProf.
                        % Otherwise, assume 'irr' was scaled. We also do this
                        Nirrs(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc(seg);
                        Nirrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc_std(seg);
                    else
                        Nirrs(Zs(j),ts(j)) = decays_fits_struct{j}.irr(seg);
                        Nirrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irr_std(seg);
                    end
                    Ntimestp(Zs(j),ts(j)) =  decays_fits_struct{j}.timestp;
                elseif strfind(zname,'FAD')
                    Fdecays{Zs(j),ts(j)} = decay;
                    Ftau1s(Zs(j),ts(j)) = restab(3,1);
                    Ftau1s_std(Zs(j),ts(j)) = restab(3,2);
                    Ftau2s(Zs(j),ts(j)) = restab(3,1)*restab(5,1);
                    Ftau2s_std(Zs(j),ts(j)) = restab(3,1)*restab(5,2);
                    Ffbound(Zs(j),ts(j)) = restab(4,1);
                    Ffbound_std(Zs(j),ts(j)) = restab(4,2);
                    ind = find(decays_fits_struct{j}.selected_pixel);
                    Fpix{Zs(j),ts(j)} = double(decays_fits_struct{j}.image(ind));
                    Fnumpix(Zs(j),ts(j)) = length(Fpix{Zs(j),ts(j)});
                    Fphot(Zs(j),ts(j)) = sum(Fpix{Zs(j),ts(j)});
                    if isfield(decays_fits_struct{j},'irrSc')
                        Firrs(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc(seg);
                        Firrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irrSc_std(seg);
                    else
                        Firrs(Zs(j),ts(j)) = decays_fits_struct{j}.irr(seg);
                        Firrs_std(Zs(j),ts(j)) = decays_fits_struct{j}.irr_std(seg);
                    end
                    Ftimestp(Zs(j),ts(j)) =  decays_fits_struct{j}.timestp;
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
        
        % Make tau1's the short lifetime
        if exist('Ntau1s')
            if size(Ntau1s,2)<length(tgrid)
                % If not enough rows/columns, pad with 0's to get the same
                % dimensions as Zgrid and tgrid
                Ntau1s(:,end+1:length(tgrid))=0; Ntau1s_std(:,end+1:length(tgrid))=0;
                Ntau2s(:,end+1:length(tgrid))=0; Ntau2s_std(:,end+1:length(tgrid))=0;
                Nfbound(:,end+1:length(tgrid))=0; Nfbound_std(:,end+1:length(tgrid))=0;
                Npix(:,end+1:length(tgrid))=num2cell(0);
                Nnumpix(:,end+1:length(tgrid))=0; Nphot(:,end+1:length(tgrid))=0;
                Nirrs(:,end+1:length(tgrid))=0; Nirrs_std(:,end+1:length(tgrid))=0;
                Ntimestp(:,end+1:length(tgrid))=0;
            end
            if size(Ntau1s,1)<length(Zgrid)
                Ntau1s(end+1:length(Zgrid),:)=0; Ntau1s_std(end+1:length(Zgrid),:)=0;
                Ntau2s(end+1:length(Zgrid),:)=0; Ntau2s_std(end+1:length(Zgrid),:)=0;
                Nfbound(end+1:length(Zgrid),:)=0; Nfbound_std(end+1:length(Zgrid),:)=0;
                Npix(end+1:length(Zgrid),:)=num2cell(0);
                Nnumpix(end+1:length(Zgrid),:)=0; Nphot(end+1:length(Zgrid),:)=0;
                Nirrs(end+1:length(Zgrid),:)=0; Nirrs_std(end+1:length(Zgrid),:)=0;
                Ntimestp(end+1:length(Zgrid),:)=0;
            end
        end
        if exist('Ftau1s')
            if size(Ftau1s,2)<length(tgrid)
                Ftau1s(:,end+1:length(tgrid))=0; Ftau1s_std(:,end+1:length(tgrid))=0;
                Ftau2s(:,end+1:length(tgrid))=0; Ftau2s_std(:,end+1:length(tgrid))=0;
                Ffbound(:,end+1:length(tgrid))=0; Ffbound_std(:,end+1:length(tgrid))=0;
                Fpix(:,end+1:length(tgrid))=num2cell(0);
                Fnumpix(:,end+1:length(tgrid))=0; Fphot(:,end+1:length(tgrid))=0;
                Firrs(:,end+1:length(tgrid))=0; Firrs_std(:,end+1:length(tgrid))=0;
                Ftimestp(:,end+1:length(tgrid))=0;
            end
            if size(Ftau1s,1)<length(Zgrid)
                Ftau1s(end+1:length(Zgrid),:)=0; Ftau1s_std(end+1:length(Zgrid),:)=0;
                Ftau2s(end+1:length(Zgrid),:)=0; Ftau2s_std(end+1:length(Zgrid),:)=0;
                Ffbound(end+1:length(Zgrid),:)=0; Ffbound_std(end+1:length(Zgrid),:)=0;
                Fpix(end+1:length(Zgrid),:)=num2cell(0);
                Fnumpix(end+1:length(Zgrid),:)=0; Fphot(end+1:length(Zgrid),:)=0;
                Firrs(end+1:length(Zgrid),:)=0; Firrs_std(end+1:length(Zgrid),:)=0;
                Ftimestp(end+1:length(Zgrid),:)=0;
            end
        end
        
        
        % Average over Z scans, keep time points
        % Use only non-zero values, where params were successfully calced
        % AVERAGES CALCULATED WITH WEIGHTING FACTORS = 1/SIGMA^2
        % http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
        
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
                Nphot_Z_aves{i}(k) = sum(Nphot(samZrange,k));
                Nnumpix_Z_aves{i}(k) = sum(Nnumpix(samZrange,k));
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
                Fphot_Z_aves{i}(k) = sum(Fphot(samZrange,k));
                Fnumpix_Z_aves{i}(k) = sum(Fnumpix(samZrange,k));
            end
        end
        
        if exist('Ntau1s')&exist('Ftau1s')
            redoxirr_Z_aves{i} = Nirr_Z_aves{i}./Firr_Z_aves{i};
            redoxirr_Z_aves_std{i} = Nirr_Z_aves{i}./Firr_Z_aves{i}.*sqrt((Nirr_Z_aves_std{i}./Nirr_Z_aves{i}).^2+(Firr_Z_aves_std{i}./Firr_Z_aves{i}).^2);
        end
        
        % Sometimes short and long lifetimes get switched, but I want tau1
        % to be the SHORT lifetime, so switch if this happens
        if exist('Ntau1s')
            if mean(Ntau2_Z_aves{i}(~isnan(Ntau2_Z_aves{i})))<mean(Ntau1_Z_aves{i}(~isnan(Ntau1_Z_aves{i})))
                tmp = Ntau2_Z_aves{i}; Ntau2_Z_aves{i} = Ntau1_Z_aves{i}; Ntau1_Z_aves{i} = tmp;
                tmp = Ntau2_Z_aves_std{i}; Ntau2_Z_aves_std{i} = Ntau1_Z_aves_std{i}; Ntau1_Z_aves_std{i} = tmp;
                %             tmp = Ntau2_t_aves{i}; Ntau2_t_aves{i} = Ntau1_t_aves{i}; Ntau1_t_aves{i} = tmp;
                %             tmp = Ntau2_t_aves_std{i}; Ntau2_t_aves_std{i} = Ntau1_t_aves_std{i}; Ntau1_t_aves_std{i} = tmp;
                Nfbound_Z_aves{i} = 1 - Nfbound_Z_aves{i}; %%Nfbound_t_aves{i} = 1 - Nfbound_t_aves{i};
            end
        end
        if exist('Ftau1s')
            if mean(Ftau2_Z_aves{i}(~isnan(Ftau2_Z_aves{i})))<mean(Ftau1_Z_aves{i}(~isnan(Ftau1_Z_aves{i})))
                tmp = Ftau2_Z_aves{i}; Ftau2_Z_aves{i} = Ftau1_Z_aves{i}; Ftau1_Z_aves{i} = tmp;
                tmp = Ftau2_Z_aves_std{i}; Ftau2_Z_aves_std{i} = Ftau1_Z_aves_std{i}; Ftau1_Z_aves_std{i} = tmp;
                %             tmp = Ftau2_t_aves{i}; Ftau2_t_aves{i} = Ftau1_t_aves{i}; Ftau1_t_aves{i} = tmp;
                %             tmp = Ftau2_t_aves_std{i}; Ftau2_t_aves_std{i} = Ftau1_t_aves_std{i}; Ftau1_t_aves_std{i} = tmp;
                Ffbound_Z_aves{i} = 1- Ffbound_Z_aves{i}; %Ffbound_t_aves{i} = 1- Ffbound_t_aves{i};
            end
        end
        
        % Use parameter matrices to populate time arrays
        if exist('Ntau1s')
            % Populate Time Arrays
            Ntimestp(Ntimestp==0)= nan;
            Ntm = mean(Ntimestp(:,samtrange),1,'omitnan');
            if ~tstmpBool Ntm = (Ntm - min(Ntm(find(Ntm))))*86400/60; end% min time units
            TimeArrs(i,seg).Ntime = Ntm;
            TimeArrs(i,seg).Nirr = Nirr_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Ntau1 = Ntau1_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Ntau2 = Ntau2_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Nbound = Nfbound_Z_aves{i}(samtrange); % NADH fbound
            TimeArrs(i,seg).Nphot = Nphot_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Nnumpix = Nnumpix_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Ntaumean = Ntau1_Z_aves{i}(samtrange).*(1-Nfbound_Z_aves{i}(samtrange))+Ntau2_Z_aves{i}(samtrange).*Nfbound_Z_aves{i}(samtrange);
            %Std time arrays to keep some estimate of error
            Ntm = mean(Ntimestp(:,samtrange),1,'omitnan');
            Ntm = (Ntm - min(Ntm(find(Ntm))))*86400/60; % min time units
            TimeArrs_std(i,seg).Ntime = Ntm;
            TimeArrs_std(i,seg).Nirr = Nirr_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Ntau1 = Ntau1_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Ntau2 = Ntau2_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Nbound = Nfbound_Z_aves_std{i}(samtrange); % NADH fbound
            % Prop of uncertainty for mean lifetimes is a little complicated.
            % Algrebra done elsewhere, taken from https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas
            % Simplify by assuming no covariance.
            tau1 = Ntau1_Z_aves{i}(samtrange); tau2 = Ntau2_Z_aves{i}(samtrange); fb = Nfbound_Z_aves{i}(samtrange);
            tau1std = Ntau1_Z_aves_std{i}(samtrange); tau2std = Ntau2_Z_aves_std{i}(samtrange); fbstd = Nfbound_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Ntaumean = sqrt(2.*(fbstd./fb).^2+(tau1std./tau1).^2+(tau2std./tau2).^2+tau2std.^2);
        end
        
        if exist('Ftau1s')
            % Populate Time Arrays
            Ftimestp(Ftimestp==0)= nan;
            Ftm = mean(Ftimestp(:,samtrange),1,'omitnan');
            if ~tstmpBool Ftm = (Ftm - min(Ftm(find(Ftm))))*86400/60;  end% min time units
            TimeArrs(i,seg).Ftime = Ftm;
            TimeArrs(i,seg).Firr = Firr_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Ftau1 = Ftau1_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Ftau2 = Ftau2_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Fbound = Ffbound_Z_aves{i}(samtrange); % FADH fbound
            TimeArrs(i,seg).Fphot = Fphot_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Fnumpix = Fnumpix_Z_aves{i}(samtrange);
            TimeArrs(i,seg).Ftaumean = Ftau1_Z_aves{i}(samtrange).*(1-Ffbound_Z_aves{i}(samtrange))+Ftau2_Z_aves{i}(samtrange).*Ffbound_Z_aves{i}(samtrange);
            %Std time arrays to keep some estimate of error
            Ftm = mean(Ftimestp(:,samtrange),1,'omitnan');
            Ftm = (Ftm - min(Ftm(find(Ftm))))*86400/60; % min time units
            TimeArrs_std(i,seg).Ftime = Ftm;
            TimeArrs_std(i,seg).Firr = Firr_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Ftau1 = Ftau1_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Ftau2 = Ftau2_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Fbound = Ffbound_Z_aves_std{i}(samtrange); % FADH fbound
            tau1 = Ftau1_Z_aves{i}(samtrange); tau2 = Ftau2_Z_aves{i}(samtrange); fb = Ffbound_Z_aves{i}(samtrange);
            tau1std = Ftau1_Z_aves_std{i}(samtrange); tau2std = Ftau2_Z_aves_std{i}(samtrange); fbstd = Ffbound_Z_aves_std{i}(samtrange);
            TimeArrs_std(i,seg).Ftaumean = sqrt(2.*(fbstd./fb).^2+(tau1std./tau1).^2+(tau2std./tau2).^2+tau2std.^2);
        end
        
        if exist('Ntau1s')&exist('Ftau1s')
            %         redoxirraves_wt(i) = Nirraves_wt(i)./Firraves_wt(i);
            %         redoxirraves_wt_std(i,seg) = Nirraves_wt(i)./Firraves_wt(i).*sqrt((Nirraves_wt_std(i)./Nirraves_wt(i)).^2+(Firraves_wt_std(i)./Firraves_wt(i)).^2);
            TimeArrs(i,seg).redox = Nirr_Z_aves{i}(samtrange)./Firr_Z_aves{i}(samtrange);
        end
        
        
        Ftm = mean(Ftimestp(:,samtrange));
        Ftm = (Ftm - min(Ftm(find(Ftm))))*86400/60; % min time units
        
        TimeArrs(i,seg).labels = labels{i};
        TimeArrs_std(i,seg).labels = labels{i};
    end % End of seg loop
end

if tstmpBool
    if exist('Label')
        save([path 'TimeArrs_tstmps_' Label '.mat'],'TimeArrs','TimeArrs_std');
    else
        save([path 'TimeArrs_tstmps.mat'],'TimeArrs','TimeArrs_std');
    end
else
    if exist('Label')
        save([path 'TimeArrs_' Label '.mat'],'TimeArrs','TimeArrs_std');
    else
        save([path 'TimeArrs.mat'],'TimeArrs','TimeArrs_std');
    end
end



