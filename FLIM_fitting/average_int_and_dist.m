function average_int_and_dist(jobind,lowoxy,highoxy, savepath)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};

load(join([analyzefolder,'/', savepath, '/',name,'_irr_decay_dist_in_rings.mat'], ''));
load(join([analyzefolder,'/', savepath, '/', name,'_fitted_in_rings.mat'], ''));

for i=1:length(flim_struct)
    if ~isempty(flim_struct(i).oxygen) 
        oxygen_level_list(i)=flim_struct(i).oxygen; 
    end
end

ind_high=find(oxygen_level_list>=lowoxy&oxygen_level_list<=highoxy);

 for ii=1:size(irr_mito_cell_dist,1)
      for jj=1:size(irr_mito_cell_dist,2)
          for kk=1:size(irr_mito_cell_dist,3)
                    
              if isempty(irr_mito_cell_dist{ii,jj,kk})
                 irr_mito_cell_dist{ii,jj,kk}=NaN;
              end  
              
              if isempty(dist_list{ii,jj,kk})
                 dist_list{ii,jj,kk}=NaN;
              end  
                    
          end
      end
 end

irr_mito_cell_dist_mat=cell2mat(irr_mito_cell_dist);
dist_mito_cell_dist_mat=cell2mat(dist_list);

for l=1:size(irr_mito_cell_dist_mat,2)
    for k=1:size(irr_mito_cell_dist_mat,3)
        
    irr_mito_cell_dist_mat_mean(l,k)=mean(irr_mito_cell_dist_mat(ind_high,l,k), 'omitnan');
    dist_mito_cell_dist_mat_mean(l,k)=mean(dist_mito_cell_dist_mat(ind_high,l,k), 'omitnan').*0.42;

    end
end

save(join([analyzefolder, '/', savepath, '/', name,'_int_fitted_in_rings.mat'], ''),'p_lm_mito_1','p_lm_mito_4','irr_mito_cell_dist_mat_mean','dist_mito_cell_dist_mat_mean');

end