clear all;
close all

addpath tracking_Kilfoil

ww = 6;
zm = 16;
mag = 20;

[file_name, pth_sdt, findex] = uigetfile2({'*.sdt';'*.tif'},'Pick a file');

if findex == 1 %sdt
    block=1;
    
    sdt = bh_readsetup([pth_sdt file_name]);
    ch = bh_getdatablock(sdt,block);
    ch = uint16(squeeze(sum(ch,1)));

elseif findex == 2 %tif
    ch = imread([pth_sdt file_name]);
    
elseif findex == 0 %cancelled
    return;
end

    ch = double(ch);
% figure
% imagesc(ch);
% axis image

featsize = 4;
masscut = 5;
Imin = 10;
field = 2;
barcc = 0.1;
barrg = 6;
barint = 10;
IdivRg = 5;

d=0;
MT=[];
Mrej=[];

M = feature2D(ch,1,featsize,masscut,Imin,field);

if isempty(M) | (M == -1)
    errordlg('Ch1: Failed to find local minimum. Try different setting')
    return
end

a = length(M(:,1));

for i=1:a
    if ((M(i,5)>barcc))
        Mrej =[Mrej ; M(i,:)];
        M (i,1:5)=0;
        %          end
        
    elseif ((M(i,4)>barrg))
        Mrej=[Mrej; M(i,:)];
        M(i,1:5)=0;
        %          end
        
    elseif ((M(i,3)<barint))
        Mrej=[Mrej; M(i,:)];
        M(i,1:5)=0;
    elseif((M(i,3)/M(i,4)<IdivRg))
        Mrej=[Mrej; M(i,:)];
        M(i,1:5)=0;
    end
end

%    Deleting the zero rows
M=M(M(:,1)~=0,:);
a = length(M(:,1));
MT(d+1:a+d, 1:5)=M(1:a,1:5);


figure;
imagesc(ch),colormap(gray);
axis image
hold on

% Making a circle the size of the feature around each feature.
h = DrawCircle(MT(:,1),MT(:,2),sqrt(MT(:,4)));

if( ~isempty(Mrej)>0 )
    plot( Mrej(:,1), Mrej(:,2), 'r.' );
end

axis image;

format short g
disp(M)
disp(['Kept : ' num2str(size(M,1))])
disp(Mrej)
disp(['Minimum Intensity : ' num2str(min(M(:,3)))])
disp(['Maximum Rg : ' num2str(max(M(:,4)))])
disp(['Maximum Eccentricity : ' num2str(max(M(:,5)))])

%%

xx = round(MT(:,1));
yy = round(MT(:,2));

xmin = xx-ww;
xmax = xx+ww;
ymin = yy-ww;
ymax = yy+ww;

c = [50,0,1];

mkdir(pth_sdt,'MeasureResult');

Sig = zeros(size(MT,1),1);
DSig = zeros(size(MT,1),1);
Nphoton = zeros(size(MT,1),1);

for i = 1:size(MT,1)
    x = xmin(i):xmax(i);
    y = ymin(i):ymax(i);

    GlobalVarTracking;
    ConstructGlobalVar(x,y,c);

    
    z = double(ch(y,x));
    
    
    %%
    [X,Y]=meshgrid(x,y);
    cmx = sum(X(:).*z(:))/sum(z(:));
    cmy = sum(Y(:).*z(:))/sum(z(:));
    cm = [cmx,cmy];
    
    pinit = [cmx,cmy,0.5,0.5,0,max(z(:))-min(z(:)),min(z(:))];
    pmax = pinit;
    pmin = pinit;
    pmin(1,1:2) = pmin(1:2)-2;
    pmax(1,1:2) = pmax(1:2)+2;
    pmin(1,3:4) = 0.1;
    pmax(1,3:4) = 4;
    pmin(1,5) = -180; pmax(1,5) = 180;
    pmin(1,6) = 0.8*pinit(1,6);
    pmax(1,6) = 1.2*pinit(1,6);
    pmin(7)= 0;
    pmax(7)= 1/2*max(z(:));
    
    dp = ones(1,7)*0.01;
    dp(4) = .01;
    dp(5) = .01;
    %%
    nonzero_z = z;
    nonzero_z(z==0)=1;
    sigz = sqrt(nonzero_z);
    weight = 1./sigz;
%     weight = ones(size(z,1),size(z,2));
    [pfit,X2,sigp,sigy,corr,Rsq,cvg_hst, converged] = lm2(@JointGaussianModel,pinit,x,y,z,weight,dp,pmin,pmax,c);
    
    zfit = JointGaussianModel(x,y,pfit,c);
    
    %%
    %pixel size in nanometer
    pixsize = 440*1000*40/mag/zm/size(ch,1);
    %pixsize = 173.9;
    %FWHM
    FWHM = pixsize*2.3548*mean([pfit(3) pfit(4)]);
    
    h = figure;
    
    subplot(1,2,1), imagesc(x,y,z)
    axis image;
    text(mean(x),min(y)-1.5,['# photons=',num2str(sum(z(:)))],'HorizontalAlignment','center');
    text(mean(x),min(y)-2.5,['pixsize=',num2str(pixsize),'nm'],'HorizontalAlignment','center');
    
    subplot(1,2,2), imagesc(x,y,zfit)
    axis image;
    
    fitstring=sprintf('%.2f,',pfit');
    stdfit = sprintf('+-%.2f,',sigp');
    
    text(mean(x),max(y)+2.5,'x,y,sigx,sigy,theta,h,bg','HorizontalAlignment','center');
    text(mean(x),max(y)+3.5,fitstring,'HorizontalAlignment','center');
    text(mean(x),max(y)+4.5,stdfit,'HorizontalAlignment','center');
    text(mean(x),min(y)-2,['# photons=',num2str(sum(zfit(:)))],'HorizontalAlignment','center');
    text(mean(x),min(y)-1,['FWHM=',num2str(FWHM), 'um'],'HorizontalAlignment','center');
    
    print(h,'-dpng',[pth_sdt 'MeasureResult/' file_name(1:end-4) '_particle' num2str(i,'%.2d') '.png'])
    
%    print(h,'-dpng','balhbalh.png')
    close(h)
    
    Sig(i) = pfit(3);
    DSig(i) = sigp(3);
    Nphoton(i) = sum(z(:));
        
end
%%
save([pth_sdt 'MeasureResult/' file_name(1:end-4) '.mat']);
disp(['Mean rad: ' num2str(mean(Sig))])
disp(['Std rad: ' num2str(mean(Sig))])
disp(std(Sig))

% 
% text(param(1,1),param(1,2)+1.8,'mean pfit','HorizontalAlignment','center');
% text(param(1,1),param(1,2)+1.6,fitstring,'HorizontalAlignment','center');
% text(param(1,1),param(1,2)+1.4,'std pfit','HorizontalAlignment','center');
% text(param(1,1),param(1,2)+1.2,stdfit,'HorizontalAlignment','center');
% text(param(1,1),param(1,2)-1.3,['mean xrms=',num2str(mean(xrms))],'HorizontalAlignment','center');
% text(param(1,1),param(1,2)-1.5,['std xrms=',num2str(std(xrms))],'HorizontalAlignment','center');