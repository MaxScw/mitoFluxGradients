%A Script that does two things:
%1) Takes raw data stored
%in s1-a1/Pos*_s1_a1.tif and information in s1-a1/s1_a1_Pos*_Masks.mat 
%to generate individual tiffs sorted by batch, position, cell, and frame number organized in separate 
%subfolders containing individual image names with format MLsets_Pos#_cell#_frame#. 
%User specifies the desired width and height of the cropped image and
%the prefix name "MLsets_" (or whatever the user wants to call it)

%2) creates a cell array fcell consisting of {fcellname, [abs(1-eigratio),totint,area,numhcpnts]}
%where fcellname has format MLsets_Pos*_cell*_frame# consisting of every each cell image in s1-a1 directory 
%Automatically saves to .mat file in the s1-a1/ directory with default name "embryo_data_features.mat" 

%set width and height of cropped image
width = 150;
height = 150;
%set prefix for name
prefix = 'MLsets_';

%Ask user for input directory s1-a1/
base = uigetdir();
cd(base);

%check for Mask files
mfiles = dir(['*_Masks.mat']);
fcell = {}; %cell array for storing image features of embryos [abs(1-eigratio),totint,area,numhcpnts]
if ~isempty(mfiles)
    for i=1:size(mfiles)  %number of Positions
        load(mfiles(i).name); %load a Position
        name = mfiles(i).name;
        under = strfind(name,'_');
        imgname = fullfile(base,'sorted_sdts', ['ML_' name(under(2)+1:under(3)) name(1:under(2)-1) '.tif']);
        %imgname = [base '/' 'sorted_sdts/' name(under(2)+1:under(3)) name(1:under(2)-1) '.tif'];
%         imgname = fullfile(base,'sorted_sdts', [name(under(3)+1:under(4)) name(1:under(3)-1) '.tif']);
%         imgname = [base '/' 'sorted_sdts/ML_' name(under(3)+1:under(4)) name(1:under(3)-1) '.tif']; %NAD/FAD ratio files
        for j=1:max(trajs(:,4))
            f = find(trajs(:,4)==j);
            if size(f,1) > 3  %screen out positions that contain fewer than 4 images
                try
                    for k=1:3:size(trajs(f,:))  %average 3 z-frames together
                    %if ~exist([base '/' 'sorted_sdts/' name(under(2)+1:under(3)) name(1:under(2)-1)])
                    %    mkdir([base '/' 'sorted_sdts/' name(under(2)+1:under(3)) name(1:under(2)-1)]);
                    %end
                    if ~exist(fullfile(base, 'sorted_sdts', [prefix name(under(2)+1:under(3)) name(1:under(2)-1)]))
                        mkdir(fullfile(base, 'sorted_sdts', [prefix name(under(2)+1:under(3)) name(1:under(2)-1)]));
                    end
                    M = zeros(size(Masks{trajs(f(k),3)}(trajs(f(k),4)).NL,2),size(Masks{trajs(f(k),3)}(trajs(f(k),4)).NL,1));
                    cropM = imcrop(M,[round(trajs(f(k),1)-width/2),round(trajs(f(k),2)-height/2),width,height]);
                    %cropM2 = cropM;
                    for l = 1:3
                        if size(f,1) >= k+l-1  %make sure that index does not exceed dimension of time frames
                            img = imread(imgname,trajs(f(k+l-1),3));
                            %img2 = imread(imgname2,trajs(f(k+l-1),3));
                            mask = Masks{trajs(f(k+l-1),3)}(trajs(f(k+l-1),4)).NL;
                            if ~isempty(img) & ~isempty(mask)
                                crop = imcrop(img,[round(trajs(f(k+l-1),1)-width/2),round(trajs(f(k+l-1),2)-height/2),width,height]);
                                %crop2 = imcrop(img2,[round(trajs(f(k+l-1),1)-width/2),round(trajs(f(k+l-1),2)-height/2),width,height]);
                                cropmask = imcrop(mask,[round(trajs(f(k+l-1),1)-width/2),round(trajs(f(k+l-1),2)-height/2),width,height]);
                                maskim = cropmask.*double(crop);
                                %maskim2 = cropmask.*double(crop2);
                                dimx = min([size(cropM,2) size(maskim,2)]);
                                dimy = min([size(cropM,1) size(maskim,1)]);
                                cropM = cropM(1:dimy,1:dimx) + maskim(1:dimy,1:dimx);
                                %cropM2 = cropM2 + maskim2;
                                M = M + mask; %combine 3 masks from different z-levels 
                            end
                        end
                    end
                    cropM = cropM/3.0; %average of 3 frames
                    %cropM2 = cropM2/3.0; %average of 3 frames
                    totint = sum(sum(cropM)); %total intensity of embryo
                    %totint2 = sum(sum(cropM2)); %total intensity of NAD/FAD ratio embryo image
                    cropM=uint8(cropM); %convert to uint8
                    %cropM2=uint8(cropM2); %convert to uint8
                    %if ~exist([base '/' 'sorted_sdts/' name(under(2)+1:under(3)) name(1:under(2)-1) '/' name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4))])
                    %    mkdir([base '/' 'sorted_sdts/' name(under(2)+1:under(3)) name(1:under(2)-1) '/' name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4))]);
                    %end
                    if ~exist(fullfile(base, 'sorted_sdts', [prefix name(under(2)+1:under(3)) name(1:under(2)-1)],[prefix name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4))]))
                        mkdir(fullfile(base, 'sorted_sdts', [prefix name(under(2)+1:under(3)) name(1:under(2)-1)],[prefix name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4))]));
                    end
                    %dirname = [base '/' 'sorted_sdts/' name(under(2)+1:under(2)) name(1:under(2)-1) '/' name(under(2)+1:under(2)) 'mask' num2str(trajs(f(k),4)) '/'];
                    %filename = [name(under(2)+1:under(2)) 'mask' num2str(trajs(f(k),4)) '_' num2str(trajs(f(k),3),'%03i') '.tif'];
                    dirname = fullfile(base,'sorted_sdts', [prefix name(under(2)+1:under(3)) name(1:under(2)-1)],[prefix name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4)) filesep]);
                    filename = [name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4)) '_' num2str(trajs(f(k),3),'%03i') '.tif'];
                    imwrite(cropM,[dirname filename],'tif','Compression','none');
                    %compute features for each cell
                    a = regionprops(M,'area'); a = a.Area;
                    p = regionprops(M,'perimeter'); p = p.Perimeter;
                    area = a; %area of cell
                    circularity = p^2/(4*pi*a);   %circularity of cell
   
                    %calculate Moment of Inertia of combined mask
                    [y,x] = find(M); ind = find(M);
                    CoM = [mean(x) mean(y)];
                    I = [[sum((y-CoM(2)).^2) ...
                        -sum((x-CoM(1)).*(y-CoM(2)))];...
                        [-sum((x-CoM(1)).*(y-CoM(2))) ...
                        sum((x-CoM(1)).^2)]];
                    Iwmass = [[sum((y-CoM(2)).^2.^M(ind)) ...
                        -sum((x-CoM(1)).*(y-CoM(2)).^M(ind))];...
                        [-sum((x-CoM(1)).*(y-CoM(2)).^M(ind)) ...
                        sum((x-CoM(1)).^2.^M(ind))]];
                    eigvals = eigs(I);
                    %eigvals = eigs(Iwmass);
                    eigratio = eigvals(1)/eigvals(2);
                    %fcellname = [prefix name(under(2)+1:under(2)) 'mask' num2str(trajs(f(k),4)) '_' num2str(trajs(f(k),3),'%03i')];
                    %calculate highcurvpnts of combined mask 
                    try
                        highcurvpnts = bwHighCurv(M,.04,.45,11);
                        numhcpnts = highcurvpnts(1);
                    catch
                        numhcpnts = 0;
                    end
                    
                    %hcpntsep = highcurvpnts(2);
                    fcellname = [prefix name(under(2)+1:under(3)) 'mask' num2str(trajs(f(k),4)) '_' num2str(trajs(f(k),3),'%03i')];
                    %fcell2 = {fcellname, [abs(1-eigratio),totint,area,totint2,numhcpnts,hcpntsep]};  %store 1-eigratio,totint,area of cell,totint2,numhcpnts,hcpntsep
                    fcell2 = {fcellname, [abs(1-eigratio),totint,area,numhcpnts]};
                    fcell = [fcell;fcell2];                
                    end
                catch
                end
            end
        end
    end
            
else
    print('No Mask files found!')
end

save(fullfile(base,'embryo_data_features.mat'),'fcell');