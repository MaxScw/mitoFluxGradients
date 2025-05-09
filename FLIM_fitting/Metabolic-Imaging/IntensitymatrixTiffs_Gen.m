function TotalInts = IntensitymatrixTiffs(path,WrVal,Ch2Bool,Bit16Bool)

% This routine uses the bh routines to read in .sdt files in a given
% folder, then converts them to intensity images and resaves them in a
% sub-folder.
% To read files into Matlab you need to read the setup parameters
% first:
%
% setupdata = readsetup(filename);
%
% In your case, the Dsdt is in the first datablock. Thus you can
% read the image Dsdt by
%
% Dsdt = bh_getdatablock(setupdata, 1)
%
% This returns a 3D matrix (t,x,y).
%
% For more complicated datafiles or to extract more information
% from the files the functions bh_blockinfo and bh_getmeasdesc
% can also be useful.
% INPUTS:
% WrVal: 0-> don't write tiffs, just calculate intensity
%        1-> write tiffs in sub folders called 'Ch1' (fluorescence),
%        'Ch2' (SHG)
%        2-> write tiffs to stacks in the same folder

% clear all;
% path = 'Z:\Lab\Tim\2016-02-05 BH Troubleshoot\2016-02-08 More tests\t1_FAD_BHonly_20srep';

%path of the Dsdt
if (isdir(path)&path(end)~='\') path = [path '\']; end
if ~exist('WrVal') WrVal = 1; end
if ~exist('Ch2Bool') Ch2Bool = 0; end
if ~exist('Bit16Bool') Bit16Bool = 0; end
Dsdt=dir([path '\*.sdt']);
sdtL = length(Dsdt);
block_1=1;
block_2=2;
%read Dsdt\

slashes = strfind(path,'\'); samlab = path(slashes(end-1)+1:end-1);

if WrVal==1
    [a1,b1] = mkdir([path '\IntTiffs']);
    if Ch2Bool [a1,b1] = mkdir([path '\Ch2Tiffs']); end;
end


sdt = bh_readsetup([path Dsdt(1).name]);
ch_1 = bh_getdatablock_v095(sdt,block_1);


% See if a tiff conversion was started before. If so, determine which
% tiffs have not been converted yet and convert those.
% ALSO, refill in the correspondence in 'nameinds'. If it errored out
% previously, nameinds may not have been saved.
frs2do = 1:sdtL;
donefrs = [];
Dtiffs = dir([path '\IntTiffs\*.tif']);
inds=cell(sdtL,1);
if size(ch_1,2)>1
    if ~isempty(Dtiffs)
        for k = 1:sdtL
            if strfind([Dtiffs(:).name],Dsdt(k).name(1:end-4))
                frs2do(frs2do==k)=[];
                donefrs = [donefrs k]; % For reloading tiffs to calculate ints
                Dsdt(k).name(1:end-4)
            end
        end
    end
end

Ch1failed  = [];
Ch2failed  = [];
IRfacts = [];

% If directory was partially converted, quickly load the previously
% converted tiffs just to calculate the total intensities
if ~isempty(donefrs)
    for i = donefrs
        sdt = bh_readsetup([path Dsdt(i).name]);
        meas = bh_getmeasdesc(sdt,1);
%         blk_info = bh_blockinfo(sdt);
        numscans(i) = double(meas.hist_fida_points);
        sdtendtimes(i) = datenum([sdt.Date ' ' sdt.Time]);
        AcqRng = GetPhotonCollectTRange_local([path Dsdt(i).name]);
        
%         asdf(i) = (AcqRng(2)-AcqRng(1))*86400;
        
        ch_1 = imread([path '\IntTiffs\' Dsdt(i).name(1:end-4) '.tif']);
        if Ch2Bool
            ch_2 = imread([path '\Ch2Tiffs\' Dsdt(i).name(1:end-4) '.tif']);
            TotalInts(i) = sum(sum(ch_1))+sum(sum(ch_2))/numscans(i);
            TotalPhots(i) = sum(sum(ch_1))+sum(sum(ch_2));
        else
            TotalInts(i) = sum(sum(ch_1))/numscans(i);
            TotalPhots(i) = sum(sum(ch_1));
        end
        
    end
end


for i = frs2do %ind:length(Dsdt)
    i
    sdt = bh_readsetup([path Dsdt(i).name]);
    meas = bh_getmeasdesc(sdt,1);
    numscans(i) = double(meas.hist_fida_points);
    sdtendtimes(i) = datenum([sdt.Date ' ' sdt.Time]);
    AcqRng = GetPhotonCollectTRange_local([path Dsdt(i).name]);

    if numscans(i)==0 numscans(i)=1; end
    try
        ch_1 = bh_getdatablock_v095(sdt,block_1);
        ch_1 = squeeze(sum(ch_1,1));
        
        if WrVal==1            
            filepath1=[path '\IntTiffs\' Dsdt(i).name(1:end-3) 'tif'];
            % imwrite(uint8(round(ch_1./max(max(ch_1)).*255)),filepath1,'tiff','Compression','none');
            if Bit16Bool
                imwrite(uint16(ch_1),filepath1,'tiff','Compression','none');
            else
                imwrite(uint8(ch_1),filepath1,'tiff','Compression','none');
            end
        elseif WrVal==2
            StkCh1(i).Dsdt=uint8(ch_1);
        end
    catch
        disp(['Im' num2str(i) 'Ch1 failed to save.'])
        Ch1failed = [Ch1failed i];
    end
    
    if Ch2Bool
        try
            ch_2 = bh_getdatablock_v095(sdt,block_2);
            ch_2 = squeeze(sum(ch_2,1));
            if WrVal==1                
                filepath2=[path '\Ch2Tiffs\' Dsdt(i).name(1:end-3) 'tif'];
                % imwrite(uint8(round(ch_2./max(max(ch_2)).*255)),filepath2,'tiff','Compression','none');
                if Bit16Bool
                    imwrite(uint16(ch_2),filepath2,'tiff','Compression','none');
                else
                    imwrite(uint8(ch_2),filepath2,'tiff','Compression','none');                    
                end
            elseif WrVal==2
                StkCh2(i).Dsdt=uint8(ch_2);
            end
            TotalInts(i) = sum(sum(ch_1))+sum(sum(ch_2))/numscans(i);
%             TotalIntsPerSec(i) = TotalInts(i)/;
            TotalPhots(i) = sum(sum(ch_1))+sum(sum(ch_2));
            Ch1Ints(i) = sum(sum(ch_1));
            Ch2Ints(i) = sum(sum(ch_2));
        catch
            disp(['Im' num2str(i) 'Ch2 failed to save.'])
            Ch2failed = [Ch2failed i];
        end
    else
        TotalInts(i) = sum(sum(ch_1))/numscans(i);
        TotalPhots(i) = sum(sum(ch_1));
    end
    
    if ~isempty(Ch1failed)
        disp(['Frames ' num2str(Ch1failed) ' failed']);
    end
    if ~isempty(Ch2failed)
        disp(['Frames ' num2str(Ch2failed) ' failed']);
    end
end

if WrVal==2
    StkWrite(StkCh1,[path 'IntenStk.tif'])
    if Ch2Bool
        StkWrite(StkCh2,[path 'IntenStk.tif'])
    end
end


if isdir(path)
    save([path '\TotalIntensities.mat'],'TotalInts','TotalPhots','numscans','IRfacts','sdtendtimes');
    sdtendtimes = (sdtendtimes-sdtendtimes(1))*24*60; % Units min
    IntP = figure;
    if exist('Ch1Ints')
        plot(Ch1Ints,'b');
        hold on; plot(Ch2Ints,'r')
    end
    [AX,H1,H2] = plotyy(sdtendtimes,TotalInts,sdtendtimes,TotalPhots);
    H2.LineStyle = '-.';
    legend('Phots/numscans','Phots','Tot','location','best')
    set(gca,'FontSize',12)
    xlabel('time (min)','fontsize',12); ylabel('Total Intensity','FontSize',14)
    title(samlab,'Interpreter','none')
    slashes = strfind(path,'\');
    saveas(gcf,[path path(slashes(end-1)+1:end-1) '_Ints.fig'])
%     close(IntP)
end;

fopen('all'); % List all open files
fclose('all'); % Close all open files



function AcqRng = GetPhotonCollectTRange_local(sdtfile)
% Uses metadata in an sdt file and extracts the exact time range during
% which photons were being collected. Units are absolute, datenum time
% (days)

sdt = bh_readsetup(sdtfile);
meas = bh_getmeasdesc(sdt,1);
numscans = double(meas.hist_fida_points);

collectend = datenum([sdt.Date ' ' sdt.Time]);
% stopt = double(meas.stop_time);
cellectduration = double(meas.fcs_end_time);
collectst = collectend - cellectduration/86400;

AcqRng = [collectst collectend];