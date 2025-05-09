function [res,txt] = xlsread3(file)
% Loads a multi-page excel spreadsheet into a 3D cell
% Author: Tim Sanchez 2018-07-10
% clear all;
% file = 'C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\2016-11-10 1-cell\ParamsAllSams.xls';


[status,sheets] = xlsfinfo(file);

for i = 1:length(sheets)
    [nums{i},txts{i},ress{i}]= xlsread(file,i);
    Dims(i,:) = size(ress{i});   
end

res = cell(max(Dims(:,1)),max(Dims(:,2)),length(sheets));

for i = 1:length(sheets)
    res(1:size(ress{i},1),1:size(ress{i},2),i) = ress{i};
    txt(1:size(txts{i},1),1:size(txts{i},2),i) = txts{i};
end

