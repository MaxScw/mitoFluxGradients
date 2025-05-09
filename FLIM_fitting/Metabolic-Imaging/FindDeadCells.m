function [DeadMask DeadInd] = FindDeadCells(Nim,Fim)
% Function to take images of embryos and identify dead cells. NADH and FAD
% images both necessary for this routine. Dead cells are bright in FAD and 
% dim in NADH. Divide FAD image by NADH. The dead cells look bright in
% these divided images. 
% Inputs: NADH and FAD tiffs
% Output: DeathMask is a binary mask showing regions that were identified
% as dead cells

% Nim = NADim; Fim = FADim;

G = fspecial('gaussian',[5 5],5);
Ngim = imfilter(Nim,G,'same');
Fgim = imfilter(Fim,G,'same');

FoverN = (Fgim+10)./(Ngim+10);

[masks DeadMask numeggs] = Masks_Kmeans_FLIMages(FoverN,-1,5,3,[2 3],[0 10^6],-1);
DeadInd = find(DeadMask);

% figure;imshow(Ngim,[])
% figure; imshow(DeadMask,[])
% figure; imshow(imdilate(DeadMask,strel('disk',2)));