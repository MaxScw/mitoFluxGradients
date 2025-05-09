% Quick script for concatenating Int tiff stacks across multiple
% acquisitions. E.g. if the embryo moved and you had to restart the
% acquisition.

clear all; close all;

pth = 'Z:\Lab\Xingbo\Media_Exchange\2018-06-29_High_rez_development_new_device\';
st = [];

for i = 1:10
    tfpth = [pth 'test' num2str(i) '\sorted_sdts\IntTiffs_Pos0_test' num2str(i)];
    D = dir([tfpth '\*.tif']);
    clear tfs;
    for j = 1:length(D)
        tfs(j) = tiffread2([tfpth '\' D(j).name]);
    end
    st = [st tfs];
end

StkWrite(st,[pth 'ConcatenatedPos0IntTiffs.tif'])