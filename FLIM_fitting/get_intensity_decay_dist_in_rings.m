function get_intensity_decay_dist_in_rings(jobind,ringnum,beginframe,stepsize,threshold, savepath)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};

sdtname = dir([pathtosdt,'/S*.sdt']);

load([analyzefolder,'/flim_structs/',name,'.mat']);
load([pathtoirf,'/','IRF_METAB_short_WD_40X_objective_pos1.mat']);
load([pathtoirf,'/illprof/','IllProfCal.mat']);


%threshold=0.7;
%beginframe=1;
%stepsize=1;

for i = beginframe:stepsize:size(flim_struct,2)-1
    for l = 1:length(flim_struct(i).cell)
              
        ind_mito_ind=find(flim_struct(i).mitoprob(flim_struct(i).cell(l).ind_embr)>threshold);
        flim_struct(i).cell(l).ind_mito=flim_struct(i).cell(l).ind_embr(ind_mito_ind);
        
    end
end

%ringnum=10;
%beginframe=1;
%stepsize=1;

for i = beginframe:stepsize:size(flim_struct,2)-1

    block_1=1;

    sdtsetup = bh_readsetup([pathtosdt,'/',sdtname(i).name]);
    sdtdata = bh_getdatablock_v095(sdtsetup,block_1);
    sdtimg = uint8(squeeze(sum(sdtdata,1)));

    %Obtain time vector for decay
    range = sdtsetup.SP_TAC_R*10^9;
    gain = double(sdtsetup.SP_TAC_G);
    resol = double(sdtsetup.SP_ADC_RE);
    dt = range/(gain*resol);
    time = (1:size(sdtdata,1))'*dt;

    for l = 1:length(flim_struct(i).cell)
        ind=[];
        dist=[];
        pixel_cell=zeros(size(sdtimg));
        pixel_cell(flim_struct(i).cell(l).ind_embr)=1;
        s = regionprops(pixel_cell,'centroid');
        aa = regionprops(pixel_cell,'MajorAxisLength');
        bb = regionprops(pixel_cell,'MinorAxisLength');
        radii = mean([aa.MajorAxisLength,bb.MinorAxisLength])./2;
        
        for k=1:length(flim_struct(i).cell(l).ind_mito)
            
        %disp(k);
        
        [row,col]=ind2sub(size(sdtimg),flim_struct(i).cell(l).ind_mito(k));
        ind(k)=flim_struct(i).cell(l).ind_mito(k);
        dist(k)=sqrt((row-s.Centroid(2)).^2+(col-s.Centroid(1)).^2);
        
        end
        
        
        [B,I]=sort(dist);
        binsize=floor(length(B)/ringnum);

        ringsize=radii./ringnum;
        
        for k=1:ringnum
            
        dist_list{i,l,k}=ringsize.*k-ringsize./2;
 
        %equal width of each ring
        pixel_mito_dist = zeros(size(sdtimg)); % dimension of 512 x 512
        decay_mito_cell_dist{i,l,k} = zeros(size(sdtdata,1),1);
       
        indnum_dist = find(dist>=ringsize.*(k-1) & dist<=ringsize.*k);
        pixel_mito_dist(flim_struct(i).cell(l).ind_mito(indnum_dist))=1;
        
        reshapedselectedpix_mito_dist = repmat(reshape(uint16(pixel_mito_dist),[1,size(pixel_mito_dist)]),[length(decay_mito_cell_dist{i,l,k}),1,1]);
        decay_mito_cell_dist{i,l,k} = decay_mito_cell_dist{i,l,k} + sum(sum(reshapedselectedpix_mito_dist.*sdtdata,2),3);
        
        
        %normalized intensity
        
        meas = bh_getmeasdesc(sdtsetup,1);
        numscans = double(meas.hist_fida_points);
        
        sdtimg_norm = double(sdtimg)./IllProfCal.*mean(IllProfCal(:));
    
        IntVals_norm_mito_dist = double(sdtimg_norm(flim_struct(i).cell(l).ind_mito(indnum_dist)));
        irr_norm_mito_dist = mean(IntVals_norm_mito_dist)/numscans;
        irr_norm_std_mito_dist = sqrt(sum(IntVals_norm_mito_dist))./(length(indnum_dist).*numscans);
    
        irr_mito_cell_dist{i,l,k}=irr_norm_mito_dist;
        irr_mito_cell_dist_err{i,l,k}=irr_norm_std_mito_dist;
        
        
        time = flim_struct(i).time;
    
        end  
    end
end

save(join([analyzefolder, '/', savepath, '/', name,'_irr_decay_dist_in_rings.mat'], ''),'-v7.3');

end