
function FLIMBatchParamPlot_wROIMovie_All(path,MaskLab)
% Performs 'FLIMBatchParamPlot_wROIMovie' on all positions in a path. Makes
% new directories for all the positions

if path(end)~='\'; path = [path '\']; end;
if ~exist('SHGBool') SHGBool = 1; end
if ~exist('AniBool') AniBool = 0; end
if ~exist('MaskLab') MaskLab = 'Masks'; end

slashes = strfind(path,'\');
Run = path(slashes (end-1)+1:end-1);
sdtpath = [path 'sorted_sdts\'];
Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    if ~strcmp(Dpos(i).name(1:3),'Pos')
        remove = [remove i];
    end
end

wrtdir = [path '\PlotswMovies\'];

[a,b] = mkdir(wrtdir);
Df = dir([path 'fits_*.mat']);

for i = 1:length(Df)
    src = [path Df(i).name];
    dest = [wrtdir Df(i).name(1:end-4) '\'];
    FLIMBatchParamPlot_wROIMovie(src, dest,MaskLab);
end
% 'PlotswMovies\'
% FLIMBatchParamPlot_wROIMovie