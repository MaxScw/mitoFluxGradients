function average_lifetime_and_dist(jobind,lowoxy,highoxy)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};


load([pathtoirf,'/','IRF_METAB_short_WD_40X_objective_pos1.mat']);

sdtsetup = bh_readsetup([pathtoirf,'/','IRF_METAB_short_WD_40X_objective_pos1.sdt']);
range = sdtsetup.SP_TAC_R*10^9;
gain = double(sdtsetup.SP_TAC_G);
resol = double(sdtsetup.SP_ADC_RE);
dt_irf = range/(gain*resol);

load([pathtoirf,'/','IRF_METAB_short_WD_40X_objective_pos1.mat']);

[val,ind]=max(irf);

irf_time=dt_irf.*ind;

load([analyzefolder,'/',name,'_irr_decay_dist_in_rings.mat']);

for i=1:length(flim_struct)
    if ~isempty(flim_struct(i).oxygen) 
        oxygen_level_list(i)=flim_struct(i).oxygen; 
    end
end

ind_high=find(oxygen_level_list>=lowoxy&oxygen_level_list<=highoxy);


for i=1:size(decay_mito_cell_dist,1)
    for l=1:size(decay_mito_cell_dist,2)
        for k=1:size(decay_mito_cell_dist,3)
        if ~isempty(decay_mito_cell_dist{i,l,k})
        lifetime_mean{i,l,k}=sum(double(decay_mito_cell_dist{i,l,k})'.*max(double(time'-irf_time),0))./sum(double(decay_mito_cell_dist{i,l,k}));
        else
        lifetime_mean{i,l,k}=NaN;    
        end
        end
    end
end

 for ii=1:size(lifetime_mean,1)
      for jj=1:size(lifetime_mean,2)
          for kk=1:size(lifetime_mean,3)
                    
              if isempty(lifetime_mean{ii,jj,kk})
                 lifetime_mean{ii,jj,kk}=NaN;
              end  
              
              if isempty(dist_list{ii,jj,kk})
                 dist_list{ii,jj,kk}=NaN;
              end  
                    
          end
      end
 end

lifetime_mean_mat=cell2mat(lifetime_mean);
dist_mito_cell_dist_mat=cell2mat(dist_list);

for l=1:size(lifetime_mean_mat,2)
    for k=1:size(lifetime_mean_mat,3)
        
    lifetime_mean_mat_mean(l,k)=nanmean(lifetime_mean_mat(ind_high,l,k));
    dist_mito_cell_dist_mat_mean(l,k)=nanmean(dist_mito_cell_dist_mat(ind_high,l,k)).*0.42;

    end
end

save([name,'_mean_lifetime_in_rings.mat'],'dist_mito_cell_dist_mat_mean','lifetime_mean_mat_mean');

end