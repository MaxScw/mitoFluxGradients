function fit_data_in_rings(jobind,coarse,zerooxyframe,lowoxy,highoxy, savepath)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};

sdtname = dir([pathtosdt,'/s*.sdt']);

load(join([analyzefolder, '/', savepath,'/',name,'_irr_decay_dist_in_rings.mat'], ''));
load([pathtoirf,'/','IRF_METAB_short_WD_40X_objective_pos1.mat']);
load([pathtoirf,'/illprof/','IllProfCal.mat']);

%% concentric ring analysis with binned frames and adjustable segmentation threshold
%fitting all parameters in original two-exp FLIM model

%use distance based binning
decay_mito_cell=decay_mito_cell_dist;


sdtsetup = bh_readsetup([pathtoirf,'/','IRF_METAB_short_WD_40X_objective_pos1.sdt']);
range = sdtsetup.SP_TAC_R*10^9;
gain = double(sdtsetup.SP_TAC_G);
resol = double(sdtsetup.SP_ADC_RE);
dt_irf = range/(gain*resol);

flimblock=1;
measdesc = bh_getmeasdesc(sdtsetup,flimblock);
T = mean([1/measdesc.min_sync_rate,1/measdesc.max_sync_rate])*10^9;

fit_start=10;
fit_end=247;

%coarse=1;

oxygen_level_list=[];

p_lm_mito_1=[];
p_lm_mito_4=[];

for i=1:length(flim_struct)
    if ~isempty(flim_struct(i).oxygen) 
        oxygen_level_list(i)=flim_struct(i).oxygen; 
    end
end

i=1;

for l=1:length(flim_struct(i).cell)
    for k=1:ringnum
        decay_mito_oxy_1{l,k}=zeros(length(time),1);
        for i=1:length(oxygen_level_list)
            if oxygen_level_list(i)>=lowoxy && oxygen_level_list(i)<=highoxy
                
                
                
               
                if size(decay_mito_cell_dist{i,l,k}) == 0
                    continue
                else
                    decay_mito_oxy_1{l,k}=decay_mito_oxy_1{l,k}+decay_mito_cell{i,l,k};
                end
            end
        end
    end
end



for l=1:length(flim_struct(i).cell)
    for k=1:ringnum
        decay_mito_oxy_4{l,k}=zeros(length(time),1);
        for i=length(oxygen_level_list)-zerooxyframe:length(oxygen_level_list)
            decay_mito_oxy_4{l,k}=decay_mito_oxy_4{l,k}+decay_mito_cell{i,l,k};
        end
    end
end


time = flim_struct(1).time;

for l=1:length(flim_struct(1).cell)
for k=1:ringnum
    
    p_lm_mito_ini = flim_struct(1).p_lm_mito; 

    p_min0 = [-200,0.7,0.05,0.01,0.02]';
    p_max0 = [200,1,5,1,1.0]';
    
    nexpo = 2; %two exponential fitting
    
    N_param = 2*nexpo+1;
    
    free=[0,1,1,1,1];
    dp = zeros(N_param,1);
    dp(2:N_param) = 0.01*free(2:N_param);
    dp(1) = 1*free(1);
    
    
    decay_mito=decay_mito_oxy_1{l,k};
    
    nonzero_decay_mito = decay_mito;
    nonzero_decay_mito(decay_mito==0)=1;
    sigy_mito = sqrt(nonzero_decay_mito);
    weight_mito = 1./sigy_mito(fit_start:fit_end);
   
    counts_mito = sum(decay_mito(fit_start:fit_end));
    
    MaxIter = -1;
    
    if counts_mito~=0   
    [p_lm_mito,Chi_sq,sigma_p_mito,sigma_y,corr,R2,cvg_hst,converged] = ...
    lm_fastest_debug(@lm_decay_model_fastest_debug,p_lm_mito_ini,time,decay_mito,weight_mito,dp,p_min0,p_max0,[nexpo,counts_mito,fit_start,fit_end],fit_start,fit_end,irf,MaxIter,dt_irf,T,coarse);
    p_lm_mito_1{l,k}=p_lm_mito;
    else
    p_lm_mito=[];   
    p_lm_mito_1{l,k}=p_lm_mito;
    end
    
    

    decay_mito=decay_mito_oxy_4{l,k};
    
    nonzero_decay_mito = decay_mito;
    nonzero_decay_mito(decay_mito==0)=1;
    sigy_mito = sqrt(nonzero_decay_mito);
    weight_mito = 1./sigy_mito(fit_start:fit_end);
   
    counts_mito = sum(decay_mito(fit_start:fit_end));
    
    MaxIter = -1;
    
    if counts_mito~=0   
    [p_lm_mito,Chi_sq,sigma_p_mito,sigma_y,corr,R2,cvg_hst,converged] = ...
    lm_fastest_debug(@lm_decay_model_fastest_debug,p_lm_mito_ini,time,decay_mito,weight_mito,dp,p_min0,p_max0,[nexpo,counts_mito,fit_start,fit_end],fit_start,fit_end,irf,MaxIter,dt_irf,T,coarse);
    p_lm_mito_4{l,k}=p_lm_mito;
    else
    p_lm_mito=[];   
    p_lm_mito_4{l,k}=p_lm_mito;
    end
    
      
end
end

save(join([analyzefolder, '/', savepath, '/', name,'_fitted_in_rings.mat'], ''),'p_lm_mito_1','p_lm_mito_4');

end