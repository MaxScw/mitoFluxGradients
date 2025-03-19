function jox_ring_kn(jobind, temperature, savepath)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};

ocr_conv_namestr = "ocr_conv_T"+string(temperature)+"C.mat"; 

load('NADH_cal.mat')
load(ocr_conv_namestr)
load(join([analyzefolder,'/',savepath, '/', name,'_int_fitted_kn_in_rings.mat'], ''));
    
    for l=1:size(p_lm_mito_1,1)
        for k=1:size(p_lm_mito_1,2)
            
            if ~isempty(p_lm_mito_1{l,k})
    
    bound_ratio_mito_nadh(l,k) = (p_lm_mito_1{l,k}(4)./(1-p_lm_mito_1{l,k}(4)));
    long_mito_nadh(l,k) = p_lm_mito_1{l,k}(3);
    short_mito_nadh(l,k) = p_lm_mito_1{l,k}(5).*p_lm_mito_1{l,k}(3);
    
    
    kn_ring(l,k)=(p_lm_mito_4{l,k}(4)./(1-p_lm_mito_4{l,k}(4)));
    
    nadhf(l,k)=irr_mito_cell_dist_mat_mean(l,k)./(c_s.*(short_mito_nadh(l,k)+long_mito_nadh(l,k).*bound_ratio_mito_nadh(l,k)));
    kox(l,k)=(bound_ratio_mito_nadh(l,k)-kn_ring(l,k));
    jox(l,k)=k_ub.*kox(l,k).*nadhf(l,k);
    
    else
              
    bound_ratio_mito_nadh(l,k) = NaN;
    long_mito_nadh(l,k) = NaN;
    short_mito_nadh(l,k) = NaN;
    
    nadhf(l,k)=NaN;
    kox(l,k)=NaN;
    jox(l,k)=NaN;
    
            end
    
        end
    end
    

save(join([analyzefolder, '/', savepath, '/', name,'_jox_in_rings_ring_kn.mat'], ''),'irr_mito_cell_dist_mat_mean','dist_mito_cell_dist_mat_mean','bound_ratio_mito_nadh','long_mito_nadh','short_mito_nadh','nadhf','kox','jox','kn_ring');

end