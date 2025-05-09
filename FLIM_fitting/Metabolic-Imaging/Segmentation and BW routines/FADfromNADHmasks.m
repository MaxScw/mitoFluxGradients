function [Fmasks] = FADfromNADHmasks(FADim,Nmasks,area_cuts,level_fact,circ_cut,drift)
% Take NADH mask and find the close corresponding FAD mask by doing
% thresh-blur-thresh method locally. This improves correspondence between
% NADH and FAD masks.
% -FADim: duh
% -Nmasks: cell of masks (same dims as FADim)
% -drift: Eggs have drifted a max of 'drift' pixels. This looks at a
%  sub-image that is drift pixels larger than the original mask on all
%  sides
Fmasks={};
if ~cellfun('isempty',Nmasks)
    for i = 1:length(Nmasks)
        clear masks areas nums bwes;
        [Ny Nx] = find(Nmasks{i});
        subindx = max([min(Nx)-drift,1]):min([max(Nx)+drift,size(Nmasks{i},2)]);
        subindy = max([min(Ny)-drift,1]):min([max(Ny)+drift,size(Nmasks{i},1)]);
        subFim = FADim(subindy,subindx);
        
        % Blur, threshold this sub-image around where NADH mask was
        K = wiener2(subFim,[5 5]);
        level = graythresh(K);
        bw1 = im2bw(K,level*level_fact);
        G = fspecial('gaussian',[15 15],10);
        gim = imfilter(double(bw1),G,'same');
        level = graythresh(gim);
        bw = im2bw(gim,level*level_fact);
        bw = imfill(bw,'holes'); % close all;imshowt(bw);figure;imshowt(im)
        
        % Do an iterative erode to separate touching eggs and discard the ones
        % that aren't in the center
        bwe = bw;
        SameNumBlobs = 0;
        eri = 1;
        while SameNumBlobs<10
            bwe = imerode(bwe,strel('disk',2));
            %         imshow(bwe); pause(0.1)
            [L nums(eri)] = bwlabel(bwe);
            if isempty(find(L))
                disp('All blobs eroded to black in FADfromNADHmasks')
                break;
            end
            
            clear cx cy
            % Find connected clusters of pixels and perimeters, then find cluster
            % that is closest in area to the NADH mask
            for j = 1:max(L(:))
                clear Lt;
                tempma = zeros(size(L,1),size(L,2));
                tempma(find(L==j)) = 1;
                Lt = tempma;
                Lt = FillEdgeHoles(Lt);
                [y,x] = find(L==j);
                cx(j,1) = mean(x); % Keep track of all CoM's
                cx(j,2) = j;
                cy(j,1) = mean(y);
                cy(j,2) = j;
                I = [[sum((y-cy(j,1)).^2) ...
                    -sum((x-cx(j,1)).*(y-cy(j,1)))];...
                    [-sum((x-cx(j,1)).*(y-cy(j,1))) ...
                    sum((y-cy(j,1)).^2)]];
                eigvals = eigs(I);
                eigratio = eigvals(1)/eigvals(2);
                circularity = abs(1-eigratio);
                %             areas(j) = size(x,1);
                %             if circularity<circ_cut & areas(j)>area_cuts(1)& areas(j)<area_cuts(2)
                %                 Lper = bwlabel(bwperim(Lt));
                %                 [yper,xper] = find(Lper);
                %                 masks{j} = Lt;
                %             else
                %                 masks{j} = [];
                %             end
            end
            
            % Find blob that is closest to the center
            CentDists = sqrt((cx(:,1)-length(subindx)/2).^2+(cy(:,1)-length(subindy)/2).^2);
            CentInd = find(CentDists==min(CentDists));
            CentBlob = zeros(size(L,1),size(L,2));
            CentBlob(find(L==CentInd)) = 1;
            bwe = CentBlob;
            bwes{eri}=bwe;
            % Keep eroding until you get 10 iterations without throwing out
            % blobs. Then you've found your big main central blob.
            if nums(eri)==1
                SameNumBlobs = SameNumBlobs + 1;
            else
                SameNumBlobs = 0;
            end;
            eri = eri + 1;
        end
        
        % Find first frame that got a stable num of eggs
        bwe = bwes{max([eri-SameNumBlobs-1,1])};
        Fmasks{i} = zeros(size(Nmasks{i},1),size(Nmasks{i},2));
        Fmasks{i}(subindy,subindx) = bwe;%bwes{eri-SameNumBlobs-1};
        %     [y,x] = find(bwperim(Fmasks{i}));
        %     imshowt(FADim); hold on; plot(x,y,'r.'); hold off
        %     disp('')
    end
end