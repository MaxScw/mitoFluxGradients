clear all;
close all

addpath tracking_Kilfoil

[file_name, pth_sdt, findex] = uigetfile2({'*.sdt';'*.tif'},'Pick a file','MultiSelect','on');

writefname = [file_name{1}(1:end-8),'.tiff'];
if findex == 1 %sdt
    block=1;
    
    for i = 1:length(file_name)
        sdt = bh_readsetup([pth_sdt file_name{i}]);
        ch = bh_getdatablock(sdt,block);
        ch = uint16(squeeze(sum(ch,1)));
        
        if i == 1
            imwrite(ch,writefname,'tiff')
        else
            imwrite(ch,writefname,'tiff','WriteMode','append')
        end
    end

elseif findex == 2 %tif
    ch = imread([pth_sdt file_name]);
    
elseif findex == 0 %cancelled
    return;
end

    ch = double(ch);
% figure
% imagesc(ch);
% axis image