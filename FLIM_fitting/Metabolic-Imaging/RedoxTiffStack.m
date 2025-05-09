
function RedoxTiffStack(path, AniBool)
% Look for a folder 'NAD_FAD_sorted' with corresponding (indexed) NAD and FAD sdt
% files in them. Could be z-stacks of NAD and FAD, or time series
% Read images in and integrate over time to get intensity images. Construct
% a tiff stack with NADH;FAD images side-by-side

if path(end)~='\' path = [path '\']; end;
Dnad = dir([path 'NAD_FAD_sorted\NAD*.sdt']);
Dfad = dir([path 'NAD_FAD_sorted\FAD*.sdt']);
sdtL = length(Dnad);
for i = 1:sdtL
    nm = Dnad(i).name;
    index(i) = str2num(nm(7:9));
end
block_1=1;
block_2=2;
Frfailed = [];
slashes = strfind(path,'\');
StkName = [UpOneDir(path) 'IntenStk_' path(slashes(end-1)+1:end-1) '.tif'];
% StkName = [path(1:end-1) '_NAD_FAD.tif'];

DStkName = dir(StkName);

% See if a tiff conversion was started before. If so, pick up where you
% left off
if ~isempty(DStkName)
    in = imfinfo(StkName);
    StkL = length(in);
    if sdtL>=StkL
        Stk = tiffread2(StkName);
        ind1 = StkL+1;
    end
else
    ind1 = 1;
    StkL = 0;
end

for i = ind1:sdtL
    nsdt = bh_readsetup([path 'NAD_FAD_sorted\' Dnad(i).name]);
    fsdt = bh_readsetup([path 'NAD_FAD_sorted\' Dfad(i).name]);
    ve = nsdt.Version;
    disp(['fr' num2str(i)]);
    % Only look at NADH image. If it's a fifo image, the corresponding FAD
    % better be, too.
    try
        nadch1 = bh_getdatablock(nsdt,block_1);
        nadch1 = squeeze(sum(nadch1,1));
        fadch1 = bh_getdatablock(fsdt,block_1);
        fadch1 = squeeze(sum(fadch1,1));
        if AniBool
            nadch2 = bh_getdatablock(nsdt,block_2);
            nadch2 = squeeze(sum(nadch2,1));
            fadch2 = bh_getdatablock(fsdt,block_2);
            fadch2 = squeeze(sum(fadch2,1));
            if str2num(ve(1))==3 % File is a fifo image
                Stk(i).data = uint8([(nadch1+nadch2) (fadch1+fadch2)]);
            else
                error('Data type not a fifo image');
            end
        else
            if str2num(ve(1))==3 
                Stk(i).data = uint8([(nadch1) (fadch1)]);
            else
                error('Data type not a fifo image');
            end
        end
    catch
       disp(['Im' num2str(i) ' failed.'])
       Frfailed = [Frfailed i];
    end
end

fopen('all'); % List all open files
fclose('all'); % Close all open files
try
    StkWrite(Stk,StkName);
catch
end








