function mito_area_and_dist(jobind,lowoxy,highoxy)

[~,~,datatoanalyze] = xlsread('datatoanalyze.xlsx');

name = datatoanalyze{jobind,1};
pathtosdt = datatoanalyze{jobind,2};
analyzefolder = datatoanalyze{jobind,3};
pathtoirf = datatoanalyze{jobind,4};
update=datatoanalyze{jobind,8};

load([analyzefolder,'/',name,'_irr_decay_dist_in_rings.mat']);

for i = beginframe:stepsize:size(flim_struct,2)-1

    for l = 1:length(flim_struct(i).cell)
        
        ind=[];
        dist=[];
        pixel_cell=zeros(512);
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
               
        indnum_dist = find(dist>=ringsize.*(k-1) & dist<=ringsize.*k);
        mito_size{i,l,k}=length(indnum_dist(:));
            
        end  
    end
end

for i=1:length(flim_struct)
    if ~isempty(flim_struct(i).oxygen) 
        oxygen_level_list(i)=flim_struct(i).oxygen; 
    end
end

ind_high=find(oxygen_level_list>=lowoxy&oxygen_level_list<=highoxy);

 for ii=1:size(mito_size,1)
      for jj=1:size(mito_size,2)
          for kk=1:size(mito_size,3)
                    
              if isempty(mito_size{ii,jj,kk})
                 mito_size{ii,jj,kk}=NaN;
              end    
                    
          end
      end
 end
 
 mito_size_mat=cell2mat(mito_size);

for l=1:size(mito_size_mat,2)
    for k=1:size(mito_size_mat,3)
        
    mito_size_mat_mean(l,k)=nanmean(mito_size_mat(ind_high,l,k));
    
    end
end

save([name,'_mito_area_in_rings.mat'], 'mito_size_mat_mean');

end