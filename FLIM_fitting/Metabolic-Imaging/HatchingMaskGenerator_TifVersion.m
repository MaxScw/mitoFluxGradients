function HatchingMaskGenerator_TifVersion(acqpath,Position)

cd(fullfile(acqpath,strcat('sorted_sdts\ProbMaps_Pos',num2str(Position),'_s1_a1_3Z')))
Listing=dir(pwd);
Listing=Listing(3:end);

Image=imread(fullfile(pwd,Listing(1).name));
 [X,Y]=meshgrid(1:1:size(Image,1),1:1:size(Image,2));
 
 
figure
imshow(Image(:,:,1))
h = drawline('SelectedColor','yellow');
m=(h.Position(2,2)-h.Position(1,2))/(h.Position(2,1)-h.Position(1,1));
HatchingMask=Y>(m*(X-h.Position(1,1)  ) + h.Position(1,2));

cd ..
mkdir(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z_TopCell'))
mkdir(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z_BottomCell'))
copyfile(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z'),strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z_Original'))
cd(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z'))

 disp('Now generating new tiffs')
for MaskIndex=1:length(Listing)

   
        MaskLoad=imread(fullfile(pwd,Listing(MaskIndex).name));
        
        
        
        
        BottomCellMask=uint8(double(MaskLoad).*(HatchingMask));
        cd ..
        cd(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z_BottomCell'))
        imwrite(BottomCellMask,Listing(MaskIndex).name)
        
        TopCellMask=uint8(double(MaskLoad).*(1-HatchingMask));
             cd ..
        cd(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z_TopCell'))
        imwrite(TopCellMask,Listing(MaskIndex).name)
        
        cd ..
        cd(strcat('ProbMaps_Pos',num2str(Position),'_s1_a1_3Z'))
        
 

end

       load(fullfile(acqpath,strcat('JointTrophMasks_Pos',num2str(Position),'.mat')));
       cd(acqpath)
save(strcat('Original_JointTrophMasks_Pos',num2str(Position),'.mat'),'Masks','frames','trajs','Coords')
         TopCellMasks=Masks;
         BottomCellMasks=Masks;
          disp('Now generating mask files')
        for MaskIndex1=1:size(Masks,1)
            for MaskIndex2=1:size(Masks,2)

                BottomCellMasks{MaskIndex1,MaskIndex2}.L=Masks{MaskIndex1,MaskIndex2}.L.*HatchingMask;
                BottomCellMasks{MaskIndex1,MaskIndex2}.Lper=Masks{MaskIndex1,MaskIndex2}.Lper.*HatchingMask;
                TopCellMasks{MaskIndex1,MaskIndex2}.L=Masks{MaskIndex1,MaskIndex2}.L.*(1-HatchingMask);
                TopCellMasks{MaskIndex1,MaskIndex2}.Lper=Masks{MaskIndex1,MaskIndex2}.Lper.*(1-HatchingMask);

            end
        end
        

Masks=TopCellMasks;
cd(acqpath)
save(strcat('TopCell_JointTrophMasks_Pos',num2str(Position),'.mat'),'Masks','frames','trajs','Coords')

Masks=BottomCellMasks;
cd(acqpath)
save(strcat('BottomCell_JointTrophMasks_Pos',num2str(Position),'.mat'),'Masks','frames','trajs','Coords')



end