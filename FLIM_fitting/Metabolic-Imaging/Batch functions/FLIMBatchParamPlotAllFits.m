function failed = FLIMBatchParamPlotAllFits(path)
% Run FLIMBatchParamPlot_v4 on each of the fits present in the path

% path = 'C:\Users\Tim\Documents\Academic - Research\Data\2014-09-11\s2m1_N2_NADHtl\';
files = subdir([path '\*fits_Pos*.mat']);
failed = {};
for i = 1:length(files)
    try
        FLIMBatchParamPlot_v4(files(i).name)
    catch
        failed = [failed;files(i).name];
    end
end