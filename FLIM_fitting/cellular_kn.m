function cellular_kn(jobind, savepath)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};

load(join([analyzefolder,'/',savepath, '/', name,'_irr_decay_dist_in_rings.mat'], ''));
load(join([analyzefolder,'/',savepath, '/', name,'_int_fitted_in_rings.mat'], ''));

for i=1:length(flim_struct)
    if ~isempty(flim_struct(i).oxygen)
        for l=1:length(flim_struct(i).cell)
        bound_ratio_mito(i,l)=flim_struct(i).cell(l).p_lm_mito(4)./(1-flim_struct(i).cell(l).p_lm_mito(4)); 
        end
    end
end

for l=1:size(bound_ratio_mito,2)
    
    kn_cell(l)=min(bound_ratio_mito(:,l));
  
end

save(join([analyzefolder, '/', savepath, '/', name,'_int_fitted_kn_in_rings.mat'], ''),'p_lm_mito_1','p_lm_mito_4','irr_mito_cell_dist_mat_mean','kn_cell','dist_mito_cell_dist_mat_mean');

end