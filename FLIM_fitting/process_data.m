function process_data(jobind)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
update=datatoanalyze{jobind,7};

sdtname = dir([pathtosdt,'/*.sdt']);

load([analyzefolder,'/',name,'_oxy.mat']);
load([analyzefolder,'/','IRF_METAB_short_WD_40X_objective_pos1.mat']);
load([analyzefolder,'/','IllProfCal.mat']);


last_frame=datatoanalyze{jobind,8};
dist_thresh=50;
step=datatoanalyze{jobind,9};
begin_frame=datatoanalyze{jobind,10};

dd=0;
ll=1;

for ii=1:length(flim_struct(1).cell)
    
    dd=dd+1;
   
j=0;
term_mito=0;
term_cyto=0;

    bound_ratio_mito_nadh=[];
    %with covariance for bound ratio err
    bound_ratio_mito_nadh_err=[];
    bound_fraction_mito_nadh=[];
    bound_fraction_mito_nadh_err=[];
    irr_mito_nadh=[];
    irr_mito_nadh_err=[];
    long_mito_nadh=[];
    long_mito_nadh_err=[];
    short_mito_nadh=[];
    short_mito_nadh_err=[];
    
    
    bound_ratio_cyto_nadh=[];
    %with covariance for bound ratio err
    bound_ratio_cyto_nadh_err=[];
    bound_fraction_cyto_nadh=[];
    bound_fraction_cyto_nadh_err=[];
    irr_cyto_nadh=[];
    irr_cyto_nadh_err=[];
    long_cyto_nadh=[];
    long_cyto_nadh_err=[];
    short_cyto_nadh=[];
    short_cyto_nadh_err=[];
    
    oxygen_level=[];
    
    mito_size=[];
    cyto_size=[];
    embr_size=[];
    

for i=begin_frame(ll):step(ll):last_frame(ll)
    %if mod(i,2)==1
    j=j+1;
    
    if i==begin_frame(ll)
    aa=ii;
    centroid_aa=flim_struct(i).cell(aa).centroid;

    bound_ratio_mito_nadh(j) = (flim_struct(i).cell(aa).p_lm_mito(4)./(1-flim_struct(i).cell(aa).p_lm_mito(4)));
    
    %with covariance for bound ratio err
    
    bound_ratio_mito_nadh_err(j) = bound_ratio_mito_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_mito(4).^2./flim_struct(i).cell(aa).p_lm_mito(4).^2+flim_struct(i).cell(aa).sigma_p_mito(4).^2./(1-flim_struct(i).cell(aa).p_lm_mito(4)).^2+2.*flim_struct(i).cell(aa).sigma_p_mito(4).^2./((1-flim_struct(i).cell(aa).p_lm_mito(4)).*flim_struct(i).cell(aa).p_lm_mito(4)));
    
    bound_fraction_mito_nadh(j) = flim_struct(i).cell(aa).p_lm_mito(4);
   
    bound_fraction_mito_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_mito(4);
    
    irr_mito_nadh(j) = flim_struct(i).cell(aa).irr_mito(3);
    
    irr_mito_nadh_err(j) = flim_struct(i).cell(aa).irr_mito(4);
    
    long_mito_nadh(j) = flim_struct(i).cell(aa).p_lm_mito(3);
   
    long_mito_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_mito(3);
    
    short_mito_nadh(j) = flim_struct(i).cell(aa).p_lm_mito(5).*flim_struct(i).cell(aa).p_lm_mito(3);
    
    short_mito_nadh_err(j) = short_mito_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_mito(5).^2./flim_struct(i).cell(aa).p_lm_mito(5).^2+flim_struct(i).cell(aa).sigma_p_mito(3).^2./flim_struct(i).cell(aa).p_lm_mito(3).^2);
   
    oxygen_level(j)=flim_struct(i).oxygen;
    
    mito_size(j)=length(flim_struct(i).cell(aa).ind_mito).*(0.43).^2;
    
    cyto_size(j)=length(flim_struct(i).cell(aa).ind_cyto).*(0.43).^2;
    
    embr_size(j)=length(flim_struct(i).cell(aa).ind_embr).*(0.43).^2;
    
    end
    
    if i>1
        
    dist=[];
        
    for bb=1:length(flim_struct(i).cell)
       
        centroid_bb=flim_struct(i).cell(bb).centroid;
        dist(bb)=sqrt((centroid_aa(1)-centroid_bb(1)).^2+(centroid_aa(2)-centroid_bb(2)).^2);
        if dist(bb)<=dist_thresh && term_mito~=1
        aa=bb;
        centroid_aa=centroid_bb;
        
        bound_ratio_mito_nadh(j) = (flim_struct(i).cell(aa).p_lm_mito(4)./(1-flim_struct(i).cell(aa).p_lm_mito(4)));
    
        %with covariance for bound ratio err
    
        bound_ratio_mito_nadh_err(j) = bound_ratio_mito_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_mito(4).^2./flim_struct(i).cell(aa).p_lm_mito(4).^2+flim_struct(i).cell(aa).sigma_p_mito(4).^2./(1-flim_struct(i).cell(aa).p_lm_mito(4)).^2+2.*flim_struct(i).cell(aa).sigma_p_mito(4).^2./((1-flim_struct(i).cell(aa).p_lm_mito(4)).*flim_struct(i).cell(aa).p_lm_mito(4)));
    
        bound_fraction_mito_nadh(j) = flim_struct(i).cell(aa).p_lm_mito(4);
   
        bound_fraction_mito_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_mito(4);
    
        irr_mito_nadh(j) = flim_struct(i).cell(aa).irr_mito(3);
    
        irr_mito_nadh_err(j) = flim_struct(i).cell(aa).irr_mito(4);
    
        long_mito_nadh(j) = flim_struct(i).cell(aa).p_lm_mito(3);
   
        long_mito_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_mito(3);
    
        short_mito_nadh(j) = flim_struct(i).cell(aa).p_lm_mito(5).*flim_struct(i).cell(aa).p_lm_mito(3);
    
        short_mito_nadh_err(j) = short_mito_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_mito(5).^2./flim_struct(i).cell(aa).p_lm_mito(5).^2+flim_struct(i).cell(aa).sigma_p_mito(3).^2./flim_struct(i).cell(aa).p_lm_mito(3).^2);
   
        oxygen_level(j)=flim_struct(i).oxygen;
        
        mito_size(j)=length(flim_struct(i).cell(aa).ind_mito).*(0.43).^2;
    
        cyto_size(j)=length(flim_struct(i).cell(aa).ind_cyto).*(0.43).^2;
    
        embr_size(j)=length(flim_struct(i).cell(aa).ind_embr).*(0.43).^2;
    
        
        end
      
    end
    
    if min(dist)>dist_thresh
       term_mito=1;  
    end
        
    end
        
    %end
 end



bound_ratio_nadh_1110 = bound_ratio_mito_nadh;

bound_fraction_nadh_1110 = bound_fraction_mito_nadh;

irr_nadh_1110 = irr_mito_nadh;

long_nadh_1110 = long_mito_nadh;

short_nadh_1110 = short_mito_nadh;

mito_size_1110 = mito_size;

cyto_size_1110 = cyto_size;

embr_size_1110 = embr_size;



bound_ratio_nadh_1110_err = bound_ratio_mito_nadh_err;

bound_fraction_nadh_1110_err = bound_fraction_mito_nadh_err;

irr_nadh_1110_err = irr_mito_nadh_err;

long_nadh_1110_err = long_mito_nadh_err;

short_nadh_1110_err = short_mito_nadh_err;




bound_ratio_nadh_mean_mito{dd} = [bound_ratio_nadh_1110];

bound_fraction_nadh_mean_mito{dd} = [bound_fraction_nadh_1110];

irr_nadh_mean_mito{dd}= [irr_nadh_1110];

long_nadh_mean_mito{dd} = [long_nadh_1110];

short_nadh_mean_mito{dd} = [short_nadh_1110];

mito_size_mean{dd} = [mito_size_1110];

cyto_size_mean{dd} = [cyto_size_1110];

embr_size_mean{dd} = [embr_size_1110];




bound_ratio_nadh_mean_mito_err{dd} = [bound_ratio_nadh_1110_err];

bound_fraction_nadh_mean_mito_err{dd} = [bound_fraction_nadh_1110_err];

irr_nadh_mean_mito_err{dd} = [irr_nadh_1110_err];

long_nadh_mean_mito_err{dd} = [long_nadh_1110_err];

short_nadh_mean_mito_err{dd} = [short_nadh_1110_err];

oxygen_level_all{dd}=oxygen_level-min(oxygen_level);



j=0;
for i=begin_frame(ll):step(ll):last_frame(ll)
    %if mod(i,2)==1
    j=j+1; 
    
    
    if i==begin_frame(ll)
    aa=ii;
    centroid_aa=flim_struct(i).cell(aa).centroid;
    
    
    bound_ratio_cyto_nadh(j) = (flim_struct(i).cell(aa).p_lm_cyto(4)./(1-flim_struct(i).cell(aa).p_lm_cyto(4)));
    
    %with covariance for bound ratio err
    
    bound_ratio_cyto_nadh_err(j) = bound_ratio_cyto_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_cyto(4).^2./flim_struct(i).cell(aa).p_lm_cyto(4).^2+flim_struct(i).cell(aa).sigma_p_cyto(4).^2./(1-flim_struct(i).cell(aa).p_lm_cyto(4)).^2+2.*flim_struct(i).cell(aa).sigma_p_cyto(4).^2./((1-flim_struct(i).cell(aa).p_lm_cyto(4)).*flim_struct(i).cell(aa).p_lm_cyto(4)));
    
    bound_fraction_cyto_nadh(j) = flim_struct(i).cell(aa).p_lm_cyto(4);
   
    bound_fraction_cyto_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_cyto(4);
    
    irr_cyto_nadh(j) = flim_struct(i).cell(aa).irr_cyto(3);
    
    irr_cyto_nadh_err(j) = flim_struct(i).cell(aa).irr_cyto(4);
    
    long_cyto_nadh(j) = flim_struct(i).cell(aa).p_lm_cyto(3);
   
    long_cyto_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_cyto(3);
    
    short_cyto_nadh(j) = flim_struct(i).cell(aa).p_lm_cyto(5).*flim_struct(i).cell(aa).p_lm_cyto(3);
    
    short_cyto_nadh_err(j) = short_cyto_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_cyto(5).^2./flim_struct(i).cell(aa).p_lm_cyto(5).^2+flim_struct(i).cell(aa).sigma_p_cyto(3).^2./flim_struct(i).cell(aa).p_lm_cyto(3).^2);
    
    end
    
    if i>1
        
    dist=[];
        
    for bb=1:length(flim_struct(i).cell)
       
        centroid_bb=flim_struct(i).cell(bb).centroid;
        dist(bb)=sqrt((centroid_aa(1)-centroid_bb(1)).^2+(centroid_aa(2)-centroid_bb(2)).^2);
        if dist(bb)<=dist_thresh && term_cyto~=1
        aa=bb;
        centroid_aa=centroid_bb;
        
        
        bound_ratio_cyto_nadh(j) = (flim_struct(i).cell(aa).p_lm_cyto(4)./(1-flim_struct(i).cell(aa).p_lm_cyto(4)));
    
        %with covariance for bound ratio err
    
        bound_ratio_cyto_nadh_err(j) = bound_ratio_cyto_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_cyto(4).^2./flim_struct(i).cell(aa).p_lm_cyto(4).^2+flim_struct(i).cell(aa).sigma_p_cyto(4).^2./(1-flim_struct(i).cell(aa).p_lm_cyto(4)).^2+2.*flim_struct(i).cell(aa).sigma_p_cyto(4).^2./((1-flim_struct(i).cell(aa).p_lm_cyto(4)).*flim_struct(i).cell(aa).p_lm_cyto(4)));
    
        bound_fraction_cyto_nadh(j) = flim_struct(i).cell(aa).p_lm_cyto(4);
   
        bound_fraction_cyto_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_cyto(4);
    
        irr_cyto_nadh(j) = flim_struct(i).cell(aa).irr_cyto(3);
    
        irr_cyto_nadh_err(j) = flim_struct(i).cell(aa).irr_cyto(4);
    
        long_cyto_nadh(j) = flim_struct(i).cell(aa).p_lm_cyto(3);
   
        long_cyto_nadh_err(j) = flim_struct(i).cell(aa).sigma_p_cyto(3);
    
        short_cyto_nadh(j) = flim_struct(i).cell(aa).p_lm_cyto(5).*flim_struct(i).cell(aa).p_lm_cyto(3);
    
        short_cyto_nadh_err(j) = short_cyto_nadh(j).*sqrt(flim_struct(i).cell(aa).sigma_p_cyto(5).^2./flim_struct(i).cell(aa).p_lm_cyto(5).^2+flim_struct(i).cell(aa).sigma_p_cyto(3).^2./flim_struct(i).cell(aa).p_lm_cyto(3).^2);
    
        end
        
    end
    
    if min(dist)>dist_thresh
       term_cyto=1;  
    end
        
    end
    
    %end
 end



bound_ratio_nadh_1110 = bound_ratio_cyto_nadh;

bound_fraction_nadh_1110 = bound_fraction_cyto_nadh;

irr_nadh_1110 = irr_cyto_nadh;

long_nadh_1110 = long_cyto_nadh;

short_nadh_1110 = short_cyto_nadh;



bound_ratio_nadh_1110_err = bound_ratio_cyto_nadh_err;

bound_fraction_nadh_1110_err = bound_fraction_cyto_nadh_err;

irr_nadh_1110_err = irr_cyto_nadh_err;

long_nadh_1110_err = long_cyto_nadh_err;

short_nadh_1110_err = short_cyto_nadh_err;



bound_ratio_nadh_mean_cyto{dd} = [bound_ratio_nadh_1110];

bound_fraction_nadh_mean_cyto{dd} = [bound_fraction_nadh_1110];

irr_nadh_mean_cyto{dd}= [irr_nadh_1110];

long_nadh_mean_cyto{dd} = [long_nadh_1110];

short_nadh_mean_cyto{dd} = [short_nadh_1110];



bound_ratio_nadh_mean_cyto_err{dd} = [bound_ratio_nadh_1110_err];

bound_fraction_nadh_mean_cyto_err{dd} = [bound_fraction_nadh_1110_err];

irr_nadh_mean_cyto_err{dd} = [irr_nadh_1110_err];

long_nadh_mean_cyto_err{dd} = [long_nadh_1110_err];

short_nadh_mean_cyto_err{dd} = [short_nadh_1110_err];

end



%%

oxygen_all_all_mito=[];
bound_ratio_all_all_mito=[];
irr_all_all_mito=[];
long_all_all_mito=[];
short_all_all_mito=[];

for i=1:length(oxygen_level_all)

    oxygen_all_all_mito=[oxygen_all_all_mito,oxygen_level_all{i}];
    bound_ratio_all_all_mito=[bound_ratio_all_all_mito,bound_ratio_nadh_mean_mito{i}];
%     irr_all_all_mito=[irr_all_all_mito,irr_nadh_mean_mito{i}./irr_nadh_mean_mito{i}(1)];
    irr_all_all_mito=[irr_all_all_mito,irr_nadh_mean_mito{i}];
    long_all_all_mito=[long_all_all_mito,long_nadh_mean_mito{i}];
    short_all_all_mito=[short_all_all_mito,short_nadh_mean_mito{i}];
    
end

%%

% plot averaged data

[N,edges,bin]=histcounts(oxygen_all_all_mito,50);
[counts,centers]=hist(oxygen_all_all_mito,50);
nadh_bound_average_mito=[];
nadh_bound_average_err_mito=[];
irr_average_mito=[];
irr_average_err_mito=[];
long_average_mito=[];
long_average_err_mito=[];
short_average_mito=[];
short_average_err_mito=[];

for i=1:length(centers)
    nadh_bound_average_mito(i)=mean(bound_ratio_all_all_mito(find(bin==i)));
    nadh_bound_average_err_mito(i)=std(bound_ratio_all_all_mito(find(bin==i)))./sqrt(length(find(bin==i)));
    
    irr_average_mito(i)=mean(irr_all_all_mito(find(bin==i)));
    irr_average_err_mito(i)=std(irr_all_all_mito(find(bin==i)))./sqrt(length(find(bin==i)));
    
    long_average_mito(i)=mean(long_all_all_mito(find(bin==i)));
    long_average_err_mito(i)=std(long_all_all_mito(find(bin==i)))./sqrt(length(find(bin==i)));
    
    short_average_mito(i)=mean(short_all_all_mito(find(bin==i)));
    short_average_err_mito(i)=std(short_all_all_mito(find(bin==i)))./sqrt(length(find(bin==i)));
end


%% OCR per cell

%load('cell_based_FLIM_data_separate_kn_no_norm_no_test1.mat');

load('NADH_conv_factor.mat');

for i=1:length(irr_nadh_mean_mito)
    kn=min(bound_ratio_nadh_mean_mito{i});
    cf_cell{i}=cf_mean.*short_nadh_mean_mito{i}./short_mean;
    nadhf_cell{i}=1000.*irr_nadh_mean_mito{i}./(cf_cell{i}.*(long_nadh_mean_mito{i}./short_nadh_mean_mito{i}.*bound_ratio_nadh_mean_mito{i}+1));
    nadhb_cell{i}=nadhf_cell{i}.*bound_ratio_nadh_mean_mito{i};
    mito_vol(i)=max(mito_size_mean{i})./max(embr_size_mean{i}).*4/3.*pi.*sqrt(max(embr_size_mean{i})./pi).^3;
    %ocr_cell{i}=(bound_ratio_nadh_mean_mito{i}-kn).*nadhf_cell{i}.*mito_vol(i);
    %ocr_cell_norm{i}=ocr_cell{i}./mean(ocr_cell{i}(1:10));
    ocr_cell{i}=(bound_ratio_nadh_mean_mito{i}-kn).*nadhf_cell{i};
    %ocr_cell_density_mean(i)=mean(ocr_cell_density{i}(1:10));
    %mito_to_embr(i)=max(mito_size_mean{i})./max(embr_size_mean{i});
    ocr_cell_norm{i}=ocr_cell{i};
end

%% plot averaged OCR data

ocr_cell_norm_all=[];

for i=1:length(ocr_cell_norm)
   
    ocr_cell_norm_all=[ocr_cell_norm_all,ocr_cell_norm{i}];
    
end


[N,edges,bin]=histcounts(oxygen_all_all_mito,50);
[counts,centers]=hist(oxygen_all_all_mito,50);
ocr_average=[];
ocr_average_err=[];

for i=1:length(centers)
    ocr_average(i)=mean(ocr_cell_norm_all(find(bin==i)));
    ocr_average_err(i)=std(ocr_cell_norm_all(find(bin==i)))./sqrt(length(find(bin==i)));
end

%%
load('ocr_conv.mat')

%% absolute ETC flux density
for i=1:length(ocr_cell_norm)
    
    ETC_flux_den{i}=ocr_cell_norm{i}.*ocr_conv;

end

%% plot averaged ETC flux data

ETC_flux_den_all=[];

for i=1:length(ocr_cell_norm)
   
    ETC_flux_den_all=[ETC_flux_den_all,ETC_flux_den{i}];
    
end


[N,edges,bin]=histcounts(oxygen_all_all_mito,50);
[counts,centers]=hist(oxygen_all_all_mito,50);
ETC_flux_den_average=[];
ETC_flux_den_average_err=[];

for i=1:length(centers)
    ETC_flux_den_average(i)=mean(ETC_flux_den_all(find(bin==i)));
    ETC_flux_den_average_err(i)=std(ETC_flux_den_all(find(bin==i)))./sqrt(length(find(bin==i)));
end
   
    %% plot average concentration

    % plot averaged data
    
    nadhf_average_all=[];
    nadhb_average_all=[];
    nadhox_average_all=[];
    nadhre_average_all=[];
     
    
   for i=1:length(oxygen_level_all)
       
       nadhox_cell{i}=nadhf_cell{i}.*min(bound_ratio_nadh_mean_mito{i});
       nadhre_cell{i}=nadhb_cell{i}-nadhf_cell{i}.*min(bound_ratio_nadh_mean_mito{i});
    
       nadhf_average_all=[nadhf_average_all,nadhf_cell{i}];
       nadhb_average_all=[nadhb_average_all,nadhb_cell{i}];
       nadhox_average_all=[nadhox_average_all,nadhox_cell{i}];
       nadhre_average_all=[nadhre_average_all,nadhre_cell{i}];
       
       
   end
    

[N,edges,bin]=histcounts(oxygen_all_all_mito,50);
[counts,centers]=hist(oxygen_all_all_mito,50);
nadhf_average=[];
nadhf_average_err=[];
nadhb_average=[];
nadhb_average_err=[];
nadhox_average=[];
nadhox_average_err=[];
nadhre_average=[];
nadhre_average_err=[];

for i=1:length(centers)
    nadhf_average(i)=mean(nadhf_average_all(find(bin==i)));
    nadhf_average_err(i)=std(nadhf_average_all(find(bin==i)))./sqrt(length(find(bin==i)));
    
    nadhb_average(i)=mean(nadhb_average_all(find(bin==i)));
    nadhb_average_err(i)=std(nadhb_average_all(find(bin==i)))./sqrt(length(find(bin==i)));
    
    nadhox_average(i)=mean(nadhox_average_all(find(bin==i)));
    nadhox_average_err(i)=std(nadhox_average_all(find(bin==i)))./sqrt(length(find(bin==i)));
    
    nadhre_average(i)=mean(nadhre_average_all(find(bin==i)));
    nadhre_average_err(i)=std(nadhre_average_all(find(bin==i)))./sqrt(length(find(bin==i)));
end

%% plot kcat

for i=1:length(oxygen_level_all)
    
    kn=min(bound_ratio_nadh_mean_mito{i});
    
    kcat_cell{i}=ocr_cell_norm{i}.*ocr_conv./(nadhf_cell{i}.*min(bound_ratio_nadh_mean_mito{i}));
 
end

%% plot average kcat

    kcat_all=[];
   
   for i=1:length(oxygen_level_all)
       
       kcat_all=[kcat_all,kcat_cell{i}];
       
   end
    

[N,edges,bin]=histcounts(oxygen_all_all_mito,50);
[counts,centers]=hist(oxygen_all_all_mito,50);
kcat_average=[];
kcat_average_err=[];


for i=1:length(centers)
    kcat_average(i)=mean(kcat_all(find(bin==i)));
    kcat_average_err(i)=std(kcat_all(find(bin==i)))./sqrt(length(find(bin==i)));
end


save([analyzefolder,'/',name,'_processed.mat'],'-v7.3');

end