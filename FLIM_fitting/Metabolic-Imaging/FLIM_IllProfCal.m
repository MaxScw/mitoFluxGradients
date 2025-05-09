function FLIM_IllProfCal(flname)
% Input flname of sdt acquisition of NADH or FAD. This smoothes the flname to 
% get an illumination proflname matrix that can be used to scale the eggs. 
% Important for calculating correct irradiance values for different eggs.
% flname = '\';

if ~exist('flname')
    [fn,pn] = uigetfile('C:\Metabolic Imaging Calibration Files\*.sdt','Select IllProfCal File');
    flname = [pn fn];
end

sdt = bh_readsetup([flname]);
field = bh_getdatablock(sdt,1);
field = squeeze(sum(field,1));
% % for multiple flnames in a folder... obsolete
% if (isdir(flname)&flname(end)~='\') flname = [flname '\']; end
% D=dir([flname '\*.sdt']); L = size(D,1);
% block_1=1;
% block_2=2;
% sdt = bh_readsetup([flname D(1).name]);
% field = zeros(sdt.SP_IMG_Y,sdt.SP_IMG_X);
% for i = 1:L
%     sdt = bh_readsetup([flname D(i).name]);
%     ch_1 = bh_getdatablock(sdt,block_1);
%     ch_1 = squeeze(sum(ch_1,1));
%     field = field + double(ch_1./L);
% end


%% Load non-uniform illumination field calibration flname
n = round(size(field,1));
G = fspecial('gaussian',[n n],50);
close all;
% smim = smooth2(field,50,50);
avefield2 = [field(:,n:-1:1) field field(:,end:-1:end-n+1)];
avefield2 = [avefield2(n:-1:1,:); avefield2; avefield2(end:-1:end-n+1,:)];
gim = imfilter(double(avefield2),G,'same');
% gim = bpassTS(double(avefield2),n/2,512);
IllProfCal = gim(n+1:end-n,n+1:end-n);
% gim = bpassTS(field,10,512);
figure;imshow(field,[]);
figure;imshow(IllProfCal,[]);%figure;imshow(smim,[]);
IllProfCalSc = uint8(round(IllProfCal./max(IllProfCal(:)).*255));
save([flname(1:end-4) '.mat'],'IllProfCal');
imwrite(uint8(IllProfCalSc),[flname(1:end-4) '.tif'],'tiff','Compression','none');
imwrite(uint8(field),[flname(1:end-4) '_raw.tif'],'tiff','Compression','none');
