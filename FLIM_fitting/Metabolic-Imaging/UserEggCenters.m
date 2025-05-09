
function UserEggCenters(path,BHzoom,posnum,fr)
% Tried some things to automatically detect boundaries in touching eggs,
% but ultimately, it was too hard to get it to recognize all eggs. Instead,
% Just make user specify egg centers. Then this routing rotates an egg mask
% to find the best overlap and uses that as the initial mask that is used
% for the active contours
% Do one pos at a time, usually only 1 position is present

% Versions:
% 2015-04-13: Add 'zoom' capability. If BH zoom other than 2X is used,
%  enter this zoom into the function and it will scale template image
% 2015-03-03: Make uManPosNum a string instead of a number
% 2015-01-08: Incorporate an active contour operation to just the first
% frame, then simply drift those masks in Make_Masks_TouchingEggs_Drift

% path = 'Z:\Tim\2015-04-08 Cembs Temp variation\s1_20p7C';
load([path '\name_indexes.mat']);
% ma = imread('C:\Users\Tim\Documents\Academic - Research\Data\Cemb_mask.tif');
% map = bwperim(ma);
egg = imread('C:\Users\Tim\Documents\Academic - Research\Data\Cemb_ave.tif');
D = dir([path '\*beadmasks.mat']);
if ~isempty(D)
    BeadBool = 1;
    load([path '\' D(1).name]);
else
    BeadBool = 0;
end

circrad = 35;

close all;
% Input variables initialization
if path(end)~='\' path = [path '\']; end;
if exist('BHzoom')
    if BHzoom~=-1 
    else
        BHzoom = 2;
    end
else
    BHzoom = 2;
end
if exist('posnum')
    if posnum~=-1 
    else
        posnum = 1;
    end
else
    posnum = 1;
end
if exist('fr')
    if fr~=-1 
    else
        fr = 1;
    end
else
    fr = 1;
end
% Adjust picture to account for different zooms. Template taken with 2X
% zoom
egg = imresize(egg,BHzoom/2);
slashes = strfind(path,'\');
Run = path(slashes(end-1)+1:end-1);
Dpos = dir([path '\sorted_sdts']); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    if length(Dpos(i).name)>5
        remove = [remove i];
    end
end
Dpos(remove)=[];
if ~isempty(Dpos)
    srtpath = [path 'sorted_sdts\'];
    for i = 1:size(Dpos,1)
        tifpaths{i} = [srtpath 'IntTiffs_' Dpos(i).name '_' Run '\'];
        Dtifs{i} = dir([tifpaths{i} '*.tif']);        % save figures with overlaid ROIs to do quick checks after batch processing
        ROIpaths{i} = [srtpath 'ROIsCheck_' Dpos(i).name '_' Run '\'];
        [a,b] = mkdir(ROIpaths{i});
    end
else
    srtpath = path;
end
close all;
imh = figure;


for posnum = 1:size(Dpos,1)
    uManPosNum = Dpos(posnum).name;
    frames = unique(sort([nameinds{strcmp(nameinds(:,2),uManPosNum)'&[nameinds{:,7}]>0,6}]));
    ts = unique(sort([nameinds{strcmp(nameinds(:,2),uManPosNum)',3}]))+1;
    Zs = unique(sort([nameinds{strcmp(nameinds(:,2),uManPosNum)',5}]))+1;
    set(gcf,'PaperPositionMode', 'auto')
    
    im2 = imread([tifpaths{posnum} Dtifs{posnum}(fr).name]);
    xdim = size(im2,2); ydim = size(im2,1);
    % Determine how many channels. If dual, use NADH channel
    if xdim==2*ydim
        xdim = xdim/2;
        im = im2(:,1:xdim);
    elseif xdim==ydim
        im = im2;
    else
        error('Channel problem')
    end
    
    if BeadBool
        bma = bmasks{1,frames(1)};
        ind = find(bma);
        im(ind) = 0;
    end
    
    imshow(im,[0 mean2(im)*4]);
    CursHand = datacursormode(gcf);
    datacursormode on;
    input('Define range over which average is to be taken with data cursors.  Then hit Enter')
    CursInf = getCursorInfo(CursHand);
%     load('C:\Users\Tim\Documents\Academic - Research\Data\curs.mat')
%     CursInf = curs;
    
    %% Circular masks
    imshow(im,[0 mean2(im)*4]); hold on;
    w = 100;
    impad = [ zeros(ydim+2*w,w) [zeros(w,xdim);im;zeros(w,xdim)] zeros(ydim+2*w,w)];
    l = graythresh(impad);
    bw = im2bw(impad,l*.8);
    singlemask = zeros(ydim,xdim);
    for i = 1:size(CursInf,2)
        Curs = CursInf(i).Position;
        
        mask = zeros(ydim,xdim);
        [X,Y] = meshgrid(1:ydim,1:xdim);
        X=X-Curs(1);Y=Y-Curs(2);
        Z = sqrt(X.^2+Y.^2);
        circmask = zeros(ydim,xdim);
        circmask(Z<circrad) = 1;
        [y,x] = find(bwperim(circmask));
        %         plot(x,y,'r.');
        %         pause(.1)
        CircMasks{i} = circmask;
        singlemask = singlemask + circmask;
        CoMxs(i) = mean(x);
    end
    
    
    %% Rotate the mask for each point
    for i = 1:size(CursInf,2)
        Curs = CursInf(i).Position;
        circpad = singlemask;
        circpad(find(CircMasks{i}))=0;
        circpad = [ zeros(ydim+2*w,w) [zeros(w,xdim);circpad;zeros(w,xdim)] zeros(ydim+2*w,w)];
        
        for j = 1:36
            rots(j) = (j-1)*10;
            eggr = double(imrotate(egg,rots(j))); % 10 deg rotation increments
            %         G = fspecial('gaussian',[10 10],8);
            %         eggr = imfilter(eggr,G,'same');
            eggr = eggr - mean(mean(eggr));
            xoff = round(Curs(1)-size(eggr,2)/2+w); yoff = round(Curs(2)-size(eggr,1)/2+w);
            xind = xoff:(xoff+size(eggr,2)-1);
            yind = yoff:(yoff+size(eggr,1)-1);
            imcmp = double(impad(yind,xind));
            circcmp = circpad(yind,xind);
            level = graythresh(eggr);
            eggrm = imerode(im2bw(eggr,level),strel('disk',7));
            ind = find(eggrm);
            pix = imcmp(ind);
            %             C = conv2(double(imcmp),double(eggr),'same');
            %             eggmults(j) = max(max(eggr.*imcmp));
            eggmults(j) = sum(sum(eggr.*imcmp));
            eggstd(j) = std(pix);
            circover(j) = sum(sum(circcmp&eggrm));
            circover(j) = circover(j) + 1;
            metr(j) = eggmults(j)*eggstd(j)/circover(j);
            %
            %                         subplot(2,3,1)
            %                         imshow(imcmp,[])
            %                         subplot(2,3,2)
            %                         imshow(eggr,[])
            %                         subplot(2,3,3); plot(eggmults,'r');
            %                         subplot(2,3,4); plot(eggstd,'g');
            %                         subplot(2,3,5); plot(circover,'b')
            %                         subplot(2,3,6); plot(metr,'b')
            % %                         pause(.4)
        end
        %         best = find(eggmults==max(eggmults));
        best = find(metr==max(metr));
        bestrots(i) = rots(best);
        eggr = double(imrotate(egg,bestrots(i))); % 10 deg rotation increments
        eggr = eggr - mean(mean(eggr));
        level = graythresh(eggr);
        eggrm = imerode(im2bw(eggr,level),strel('disk',15));
        xoff = round(Curs(1)-size(eggr,2)/2+w); yoff = round(Curs(2)-size(eggr,1)/2+w);
        xind = xoff:(xoff+size(eggr,2)-1);
        yind = yoff:(yoff+size(eggr,1)-1);
        % padded mask, insert, then trim
        mask = zeros(ydim+2*w,xdim+2*w);
        mask(yind,xind) = eggrm;
        mask = mask(w+1:w+ydim,w+1:w+xdim);
        [y,x] = find(bwperim(mask));
        h(i) = plot(x,y,'r.');
        text(mean(x),mean(y),num2str(i),'fontsize',13,'color','r')
        IniMasks{i} = mask;
        CoMxs(i) = mean(x);
    end
    ex = 0;
    %%
    while ex==0
        egstr = input('Rotate a mask? (egg num or ''n'')');
        if isstr(egstr)
            ex = 1;
        elseif isnumeric(egstr)
            eg = egstr;
            if eg<=size(CursInf,2)
                %             eg = 12;
                rot = input('Degrees?');
                %             rot = 20;
                bestrots(eg) = bestrots(eg) + rot;
                eggr = double(imrotate(egg,bestrots(eg))); % 10 deg rotation increments
                eggr = eggr - mean(mean(eggr));
                level = graythresh(eggr);
                eggrm = imerode(im2bw(eggr,level),strel('disk',15));
                xoff = round(CursInf(eg).Position(1)-size(eggr,2)/2+w); yoff = round(CursInf(eg).Position(2)-size(eggr,1)/2+w);
                xind = xoff:(xoff+size(eggr,2)-1);
                yind = yoff:(yoff+size(eggr,1)-1);
                % padded mask, insert, then trim
                mask = zeros(ydim+2*w,xdim+2*w);
                mask(yind,xind) = eggrm;
                mask = mask(w+1:w+ydim,w+1:w+xdim);
                [y,x] = find(bwperim(mask));
                IniMasks{eg} = mask;
                CoMxs(eg) = mean(x);
                imshow(im,[0 mean2(im)*4]); hold on;
                for j = 1:length(IniMasks)
                    [y,x] = find(bwperim(IniMasks{j}));
                    plot(x,y,'r.')
                    text(mean(x),mean(y),num2str(j),'fontsize',13,'color','r')
                end
            else
                disp('There aren''t that many eggs, dumbass')
            end
        end
    end
    
    % ACTIVE CONTOURS
    disp('BEHOLD, active contours')
    [acs acscoords] = EllMasks_TouchingEggs_ActiveConts(im,IniMasks,75,6);
    ex = 0;
    imshow(im,[0 mean2(im)*4]); hold on;
    for j = 1:length(acs)
        [y,x] = find(bwperim(acs{j}));
        plot(x,y,'r.')
        text(mean(x),mean(y),num2str(j),'fontsize',13,'color','r')
    end
    while ex==0
        egstr = input('Adjust an active contour? (egg num or ''n'')');
        if isstr(egstr)
            ex = 1;
        elseif isnumeric(egstr)
            eg = egstr;
            if eg<=size(CursInf,2)
                niter = input('Iterations?');
                sm = input('Smooth Factor?');
                [acs2 acscoords] = EllMasks_TouchingEggs_ActiveConts(im,IniMasks,niter,sm,eg);
                acs{eg} = acs2{eg};
                % Remove overlap with previous masks:
                over = zeros(ydim,xdim);
                for j = 1:length(acs)    
                    if j ~= eg
                        over = over|(acs{eg}&acs{j});
                    end
                end
                acs{eg} = acs{eg} - over;
                imshow(im,[0 mean2(im)*4]); hold on;
                for j = 1:length(acs)
                    [y,x] = find(bwperim(acs{j}));
                    plot(x,y,'r.')
                    text(mean(x),mean(y),num2str(j),'fontsize',13,'color','r')
                end
            else
                disp('There aren''t that many eggs, dumbass')
            end
        end
    end
    
%     [a,b] = sort(CoMxs);
%     IniMasks = IniMasks(b);
    IniMasks = acs;
    save([path path(slashes(end-1)+1:end-1) '_' uManPosNum '_IniMasks.mat'],'frames','IniMasks');
end